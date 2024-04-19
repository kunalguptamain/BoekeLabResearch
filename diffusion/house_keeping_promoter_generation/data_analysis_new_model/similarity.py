import json
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from multiprocessing import Pool, cpu_count
import warnings

warnings.filterwarnings("ignore")

def calculate_alignment(args):
    generated_sequence, check = args
    alignments = pairwise2.align.localxx(generated_sequence, check)
    best_alignment = max(alignments, key=lambda x: x[2])  # Get the alignment with the best score
    return best_alignment

if __name__ == '__main__':
    generated_sequences_path = r'C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_new_model/new_promoters.txt'
    all_sequences_path = r'C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data/compare_sequences.txt'

    with open(generated_sequences_path, 'r') as file:
        generated_sequences = [line.strip() for line in file.readlines()]

    with open(all_sequences_path, 'r') as file:
        all_sequences = json.load(file)
        human_house_keeping_promoters = [promoter for promoter, label in all_sequences if label == 0]

    print(len(generated_sequences))
    print(len(human_house_keeping_promoters))

    generated_sequence = generated_sequences[6]

    # Create a pool of workers
    with Pool(cpu_count()) as p:
        # Calculate alignments in parallel
        resultsA = p.map(calculate_alignment, [(generated_sequence, check) for check in human_house_keeping_promoters])
    with Pool(cpu_count()) as p:
        # Calculate alignments in parallel
        resultsB = p.map(calculate_alignment, [(str(Seq(generated_sequence).reverse_complement()), check) for check in human_house_keeping_promoters])
    results = resultsA = resultsB
    # Find the best alignment
    best_alignment = max(results, key=lambda x: x[2])  # Get the alignment with the best score

    s1_list = list(best_alignment[0])
    s2_list = list(best_alignment[1])
    print(best_alignment[3])

    # Replace characters in s1 with dashes where s2 has dashes
    for i in range(len(s1_list)):
        if s2_list[i] == "-":
            s1_list[i] = "-"

    # Replace characters in s2 with dashes where s1 has dashes
    for i in range(len(s2_list)):
        if s1_list[i] == "-":
            s2_list[i] = "-"

    # Convert lists back to strings
    s1_dash = "".join(s1_list)
    print(s1_dash)
