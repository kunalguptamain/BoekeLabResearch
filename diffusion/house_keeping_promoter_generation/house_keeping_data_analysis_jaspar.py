from pyjaspar import jaspardb
import numpy as np
import csv
from Bio.Seq import Seq
import json
from enum import Enum

threshold_dissimilarity = 0.1

promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data/compare_sequences.txt"
factor_txt_dict_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/factors_{threshold_dissimilarity * 100}_dict.txt"
factor_txt_total_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/factors_{threshold_dissimilarity * 100}_total.txt"
factor_csv_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/factors_{threshold_dissimilarity * 100}.csv"

jdb_obj = jaspardb(release='JASPAR2024')

motifs = jdb_obj.fetch_motifs(
    collection = ['CORE'],
    all_versions = False,
    tax_group='Vertebrates',
)

one_hot_encode_mapping = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1],
}

def count_single_motif(sequence, target_motif, threshold):
    pssm = target_motif.pwm.log_odds()
    search = list(pssm.search(Seq(sequence), threshold=3.0))
    return len(search)

def get_motifs(sequence, threshold_dissimilarity, class_dict, general_dict):
    for i, motif in enumerate(motifs):
        count = count_single_motif(sequence, motif, threshold_dissimilarity)
        if(count == 0): continue
        if motif.name not in class_dict: class_dict[motif.name] = 0
        class_dict[motif.name] += count
        if motif.name not in general_dict: general_dict[motif.name] = 0
        general_dict[motif.name] += count

def one_hot_encode(seq):
    encoded = [one_hot_encode_mapping[i] for i in seq]
    return np.array(encoded).T

dicts = [{}, {}, {}, {}, {}, {}, {}]
total = [0, 0, 0, 0, 0, 0, 0]

def save(object, path):
    with open(path, 'w') as file:
        file.write(json.dumps(object))

with open(promoter_abs_path) as file:
    sequence_set = json.load(file)

i = 0
for sequence, label in sequence_set:
    print(i)
    get_motifs(sequence, threshold_dissimilarity, \
        dicts[label], dicts[-1])
    total[label] += 1
    total[-1] += 1
    i += 1

save(dicts, factor_txt_dict_save_path)
save(total, factor_txt_total_save_path)
# with open(factor_txt_dict_save_path, 'r') as file:
#     dicts = json.load(file)
# with open(factor_txt_total_save_path, 'r') as file:
#     total = json.load(file)

print(total)
factors = list(dicts[-1].keys())
counts = list(dicts[-1].values())

sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i])

sorted_factors = [factors[i] for i in sorted_indices]
sorted_counts = [counts[i] for i in sorted_indices]

label_dict = {
    0: "HUMAN_CODING",
    1: "HUMAN_NON_CODING",
    2: "MOUSE_CODING",
    3: "MOUSE_NON_CODING",
    4: "NON_PROMOTER",
    5: "HUMAN_CODING_GENERATED",
    6: "GENERAL"
}

with open(factor_csv_save_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([''] + [f'{label_dict[i]}' for i in range(len(dicts))])
    for key in sorted_factors:
        row = [key] + [d.get(key, 0) / total[i] for i, d in enumerate(dicts)]
        writer.writerow(row)