from standardized_dna_diffusion_pipeline.sequence_set import SequenceSet
from Bio.Seq import Seq
from pydantic import BaseModel
from typing import List, Dict, Callable
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from matplotlib.colors import to_rgba
from functools import partial
import multiprocessing

class Analysis(BaseModel):
    """
    Pydantic model for a single analysis.
    """
    name: str
    positions: List[int]
    magnitudes: List[float]


class AnalysisPerSequence(BaseModel):
    """
    Pydantic model that maps a sequence to a list of analyses.
    """
    sequence_id: int
    analyses: List[Analysis] = []

    def add_analysis(self, analysis: Analysis):
        self.analyses.append(analysis)


class Analyzer:
    def __init__(self, sequence_set: SequenceSet):
        self.sequence_set = sequence_set
        self.analysis_results: Dict[int, AnalysisPerSequence] = {}

    def query(self, condition: Callable[[SequenceSet.Sequence], bool]) -> List[SequenceSet.Sequence]:
        """
        General query function that takes a condition (lambda) and returns matching sequences.
        """
        return [sequence for sequence in self.sequence_set.sequences if condition(sequence)]

    def query_by_class(self, class_type: str) -> List[SequenceSet.Sequence]:
        """
        Queries sequences by their class type.
        """
        return self.query(lambda seq: seq.dna_type == class_type)

    def query_by_id(self, sequence_id: int) -> List[SequenceSet.Sequence]:
        """
        Queries sequences by their ID.
        """
        return self.query(lambda seq: seq.id == sequence_id)

    def query_all(self) -> List[SequenceSet.Sequence]:
        """
        Returns all sequences.
        """
        return self.query(lambda seq: True)  # Return all sequences by using a condition that is always True

    def analyze_gc_content(self, sequence: SequenceSet.Sequence) -> Analysis:
        """
        Analysis function for GC content.
        Takes in a single sequence and returns positions and magnitudes.
        """
        gc_content_positions = [i for i, base in enumerate(sequence.sequence) if base in ['G', 'C']]
        magnitudes = [1.0] * len(gc_content_positions)  # Default magnitude is 1
        return Analysis(name="GC Content", positions=gc_content_positions, magnitudes=magnitudes)

    def analysis_tfs(self, sequence: SequenceSet.Sequence, motifs, threshold: float = 0.8, name = "TF Binding Site") -> Analysis:
        """
        General transcription factor binding site analysis using motifs from the JASPAR database.
        Takes in a sequence and a list of motifs, and searches for binding sites using PWM with a given threshold.
        
        :param sequence: The DNA sequence to analyze.
        :param motifs: A list of motifs (from JASPAR or similar) to search for.
        :param threshold: Minimum PWM log-odds score to consider a match (default = 0.8).
        :return: Analysis object containing positions and magnitudes for TF binding sites.
        """
        tf_binding_positions = []
        magnitudes = []

        # Loop over each motif and search for binding sites in the sequence
        for motif in motifs:
            pwm = motif.pssm  # Position-Specific Scoring Matrix (PSSM)
            
            # Search for motif occurrences in the sequence using the given threshold
            matches = pwm.search(Seq(sequence.sequence), threshold=threshold)

            for match in matches:
                position = match[0]  # Position of the motif occurrence
                score = match[1]     # Score of the match (log-odds)
                
                # Store position and score
                tf_binding_positions.append(position)
                magnitudes.append(score)

        # Return analysis result with positions and magnitudes
        return Analysis(name=name, positions=tf_binding_positions, magnitudes=magnitudes)

    def perform_analysis(
        self, 
        analysis_func: Callable[[SequenceSet.Sequence], Analysis], 
        sequences: List[SequenceSet.Sequence], 
        params: Dict = {},
    ):
        """
        General function that performs an analysis on a list of sequences using multithreading
        and updates the results.
        
        :param analysis_func: The analysis function to apply to each sequence.
        :param sequences: List of SequenceSet.Sequence objects to analyze.
        :param max_workers: Number of threads to use (default=4).
        """
        max_workers = multiprocessing.cpu_count()
        
        def analyze_sequence(seq: SequenceSet.Sequence):
            analysis_result = analysis_func(seq, **params)

            # Thread-safe way to update the shared data structure
            if seq.id not in self.analysis_results:
                self.analysis_results[seq.id] = AnalysisPerSequence(sequence_id=seq.id)

            self.analysis_results[seq.id].add_analysis(analysis_result)

        # Use ThreadPoolExecutor to parallelize the analysis
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(analyze_sequence, sequences)

    def plot_analysis(self, queries: List[List[SequenceSet.Sequence]], attributes: List[str], query_labels: List[str], normalize: bool = True):
        """
        Plots the sum of magnitudes for selected attributes, normalized per attribute,
        so that the total height per attribute equals 1 for each query.

        :param queries: List of lists of sequences, each representing a separate query result.
        :param attributes: List of attributes (e.g., "GC Content", "TF Binding Site").
        :param query_labels: Labels for each query, used in the plot legend.
        :param normalize: If True, normalize magnitudes so that the total height = 1; if False, divide by the total length of data.
        """
        total_positions = self.sequence_set.sequence_padded_length 

        if len(queries) != len(query_labels):
            raise ValueError("The number of query groups must match the number of query labels.")

        # Prepare a base color for each query using a colormap
        base_colors = plt.cm.viridis(np.linspace(0, 1, len(queries)))

        # Loop over each query group (list of sequences)
        for query_idx, sequences in enumerate(queries):
            # Initialize a dictionary to store cumulative magnitudes for each attribute in this query group
            attribute_magnitude_sums = {attr: np.zeros(total_positions) for attr in attributes}

            for seq in sequences:
                # Get the analyses for this sequence
                if seq.id not in self.analysis_results:
                    continue  # Skip sequences without any analyses

                analysis_per_sequence = self.analysis_results[seq.id]

                # Loop through each analysis and process if it matches the requested attributes
                for analysis in analysis_per_sequence.analyses:
                    if analysis.name in attributes:
                        for i, pos in enumerate(analysis.positions):
                            attribute_magnitude_sums[analysis.name][pos] += analysis.magnitudes[i]

            # Normalize each attribute's magnitudes or divide by the total length based on the `normalize` flag
            for attr in attributes:
                total_magnitude = np.sum(attribute_magnitude_sums[attr])
                if normalize:
                    if total_magnitude > 0:
                        attribute_magnitude_sums[attr] /= total_magnitude  # Normalize to sum to 1
                else:
                    attribute_magnitude_sums[attr] /= float(len(sequences))  # Average by total sequences

            # Plot each attribute's magnitudes for this query group
            for attr_idx, (attr, magnitudes) in enumerate(attribute_magnitude_sums.items()):
                # Slightly modify the base color for different attributes
                color = np.array(to_rgba(base_colors[query_idx]))  # Get the base color for the query
                color[0] = min(1, color[0] + 0.1 * attr_idx)  # Adjust the red channel slightly
                color[1] = min(1, color[1] + 0.1 * attr_idx)  # Adjust the green channel slightly
                color[2] = min(1, color[2] + 0.1 * attr_idx)  # Adjust the blue channel slightly
                
                # Plot the results
                plt.plot(
                    range(total_positions),
                    magnitudes,
                    label=f"{query_labels[query_idx]} - {attr} ({'normalized' if normalize else 'averaged'})",
                    color=color
                )

        plt.xlabel("Position")
        plt.ylabel("Magnitude")
        plt.legend()
        plt.show()
