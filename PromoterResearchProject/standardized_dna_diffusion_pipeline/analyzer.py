from standardized_dna_diffusion_pipeline.sequence_set import SequenceSet
from Bio.Seq import Seq
from pydantic import BaseModel
from typing import List, Dict, Callable
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from matplotlib.colors import to_rgba
from Bio import motifs
import multiprocessing
from tqdm import tqdm
from collections import defaultdict

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

class FrequentMotifs:
    def __init__(self):
        pass

    def create_count_matrix(sequence):
        counts = defaultdict(lambda: [0] * len(sequence))  # Initialize counts
        for i, nucleotide in enumerate(sequence):
            if nucleotide == 'A':
                counts['A'][i] += 1
            elif nucleotide == 'C':
                counts['C'][i] += 1
            elif nucleotide == 'G':
                counts['G'][i] += 1
            elif nucleotide == 'T':
                counts['T'][i] += 1
            elif nucleotide == 'W':
                counts['A'][i] += 0.5
                counts['T'][i] += 0.5
            elif nucleotide == 'R':
                counts['A'][i] += 0.5
                counts['G'][i] += 0.5
            elif nucleotide == 'S':
                counts['C'][i] += 0.5
                counts['G'][i] += 0.5
            elif nucleotide == 'Y':
                counts['C'][i] += 0.5
                counts['T'][i] += 0.5
            elif nucleotide == 'K':
                counts['G'][i] += 0.5
                counts['T'][i] += 0.5
            elif nucleotide == 'M':
                counts['A'][i] += 0.5
                counts['C'][i] += 0.5
            elif nucleotide == 'B':
                counts['C'][i] += 0.33
                counts['G'][i] += 0.33
                counts['T'][i] += 0.33
                
        return counts

    initiator_motif = motifs.Motif(counts={
        'A': [49, 0, 288, 26, 77, 67, 45, 50],
        'C': [48, 303, 0, 81, 95, 118, 85, 96],
        'G': [69, 0, 0, 116, 0, 46, 73, 56],
        'T': [137, 0, 15, 80, 131, 72, 100, 101]
    })
    tata_motif = motifs.Motif(counts={
        'A': [61, 16, 352, 3, 354, 268, 360, 222, 155, 56, 83, 82, 82, 68, 77],
        'C': [145, 46, 0, 10, 0, 0, 3, 2, 44, 135, 147, 127, 118, 107, 101],
        'G': [152, 18, 2, 2, 5, 0, 20, 44, 157, 150, 128, 128, 128, 139, 140],
        'T': [31, 309, 35, 374, 30, 121, 6, 121, 33, 48, 31, 52, 61, 75, 71]
    })
    ccaat_motif = motifs.Motif(counts={
        'A': [56, 32, 25, 102, 51, 0, 0, 175, 119, 17, 23, 116],
        'C': [55, 52, 47, 1, 6, 173, 174, 0, 8, 0, 90, 6],
        'G': [12, 43, 24, 70, 99, 1, 0, 0, 21, 15, 59, 52],
        'T': [52, 48, 79, 2, 19, 1, 1, 0, 27, 143, 3, 1]
    })
    gc_motif = motifs.Motif(counts={
        'A': [102, 97, 50, 67, 0, 2, 54, 46, 1, 79, 23, 0, 20, 40],
        'C': [40, 31, 6, 1, 0, 0, 170, 1, 3, 0, 17, 166, 86, 24],
        'G': [50, 112, 154, 206, 274, 272, 0, 224, 222, 171, 192, 35, 52, 109],
        'T': [82, 34, 64, 0, 0, 0, 50, 3, 48, 24, 42, 73, 116, 101]
    })
    rbm =  motifs.Motif(counts=create_count_matrix("CTGGGARWTGTAGTY"))
    dre =  motifs.Motif(counts=create_count_matrix("WATCCGATW"))
    breu =  motifs.Motif(counts=create_count_matrix("SSRCGCC"))


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

    def analysis_tfs(self, sequence: SequenceSet.Sequence, motifs: List, input_name: str, threshold: float = 0.8) -> Analysis:
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
        return Analysis(name=input_name, positions=tf_binding_positions, magnitudes=magnitudes)

    def perform_analysis(
        self, 
        analysis_func: Callable[[SequenceSet.Sequence], Analysis], 
        sequences: List[SequenceSet.Sequence],
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
            analysis_result = analysis_func(seq)
            # Thread-safe way to update the shared data structure
            if seq.id not in self.analysis_results:
                self.analysis_results[seq.id] = AnalysisPerSequence(sequence_id=seq.id)

            self.analysis_results[seq.id].add_analysis(analysis_result)

        # Use ThreadPoolExecutor to parallelize the analysis
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            list(tqdm(executor.map(analyze_sequence, sequences), total=len(sequences)))

    def plot_analysis(self, queries: List[List[SequenceSet.Sequence]], attributes: List[str], query_labels: List[str], normalize: bool = True, reduce_noise: bool = False, smoothing_window: int = 4):
        """
        Plots the sum of magnitudes for selected attributes, normalized per attribute,
        so that the total height per attribute equals 1 for each query.

        :param queries: List of lists of sequences, each representing a separate query result.
        :param attributes: List of attributes (e.g., "GC Content", "TF Binding Site").
        :param query_labels: Labels for each query, used in the plot legend.
        :param normalize: If True, normalize magnitudes so that the total height = 1; if False, divide by the total length of data.
        :param reduce_noise: If True, apply a basic noise reduction to smooth out fluctuations.
        :param smoothing_window: Window size for noise reduction. Larger values smooth more.
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

                # Apply noise reduction (e.g., a moving average filter)
                if reduce_noise:
                    attribute_magnitude_sums[attr] = self.smooth_series(attribute_magnitude_sums[attr], window=smoothing_window)

            # Plot each attribute's magnitudes for this query group
            for attr_idx, (attr, magnitudes) in enumerate(attribute_magnitude_sums.items()):
                # Slightly modify the base color for different attributes
                color = np.array(to_rgba(base_colors[query_idx]))  # Get the base color for the query
                color[0] = min(1, color[0] + 0.2 * attr_idx)  # Adjust the red channel slightly
                color[1] = min(1, color[1] + 0.2 * attr_idx)  # Adjust the green channel slightly
                color[2] = min(1, color[2] + 0.2 * attr_idx)  # Adjust the blue channel slightly
                
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
    def smooth_series(self, data: np.ndarray, window: int = 5) -> np.ndarray:
        """
        Applies a moving average to smooth the data, preserving the min and max values.

        :param data: Array of magnitudes.
        :param window: Size of the moving average window.
        :return: Smoothed data array with min and max values preserved.
        """
        if len(data) < window:
            return data  # If the data length is smaller than the window, no smoothing

        smoothed_data = np.copy(data)
        half_window = window // 2

        for i in range(half_window, len(data) - half_window):
            # Preserve the original min and max
            min_val = np.min(data)
            max_val = np.max(data)
            # Apply moving average within the window
            smoothed_data[i] = np.mean(data[i - half_window:i + half_window + 1])

            # Reassign original min and max
            smoothed_data = np.clip(smoothed_data, min_val, max_val)

        return smoothed_data