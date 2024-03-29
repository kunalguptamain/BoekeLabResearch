from pyjaspar import jaspardb
import numpy as np
import csv
from Bio.Seq import Seq
import json

threshold_dissimilarity = 0.1

promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/strong_weak_promoters/promoter.txt"
non_promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/strong_weak_promoters/non_promoter.txt"
generated_promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/strong_weak_promoters/generated_promoters.txt"
factor_txt_dict_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/data_analysis_jaspar/factors_{threshold_dissimilarity * 100}_dict.txt"
factor_txt_total_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/data_analysis_jaspar/factors_{threshold_dissimilarity * 100}_total.txt"
factor_csv_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/data_analysis_jaspar/factors_{threshold_dissimilarity * 100}.csv"
jdb_obj = jaspardb(release='JASPAR2024')

motifs = jdb_obj.fetch_motifs(
    collection = ['CORE'],
    all_versions = False,
    tax_group='Escherichia coli K-12'
)

cycle = 0
save_threshold = 30
save_cycle = 0

classification_dict = {
    "ALL": 4,
    "GENERATED": 3,
    "STRONG": 2,
    "WEAK": 1,
    "NON": 0,
}

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
    print(cycle)
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

dicts = [{}, {}, {}, {}, {}]
total = [0, 0, 0, 0, 0]

if(save_cycle != 0):
    dicts = open(factor_txt_dict_save_path)
    total = open(factor_txt_total_save_path)

def save(object, path):
    with open(path, 'w') as file:
        file.write(json.dumps(object))

with open(promoter_abs_path, 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        if(i % 2 == 0):
            placed = True
            tokens = line.strip().split()
            forward = tokens[2] == "FORWARD"
            if(tokens[-1] not in classification_dict):
                placed = False
                continue
            classification = classification_dict[tokens[-1]]
        else:
            cycle += 1
            if(cycle < save_cycle): continue
            if(not placed): continue
            sequence = line.strip()
            get_motifs(sequence, threshold_dissimilarity, \
                dicts[classification], dicts[-1])
            total[classification] += 1
            total[-1] += 1
            if(cycle % save_threshold == 0):
                save(dicts, factor_txt_dict_save_path)
                save(total, factor_txt_total_save_path)
            

with open(non_promoter_abs_path, 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        cycle += 1
        if(cycle < save_cycle): continue
        sequence = line.strip()
        get_motifs(one_hot_encode(sequence), threshold_dissimilarity, \
                dicts[classification], dicts[-1])
        total[classification] += 1
        total[-1] += 1
        if(cycle % save_threshold == 0):
            save(dicts, factor_txt_dict_save_path)
            save(total, factor_txt_total_save_path)

with open(generated_promoter_abs_path, 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        cycle += 1
        if(cycle < save_cycle): continue
        sequence = line.strip()
        get_motifs(one_hot_encode(sequence), threshold_dissimilarity, \
                dicts[classification], dicts[-1])
        total[classification] += 1
        total[-1] += 1
        if(cycle % save_threshold == 0):
            save(dicts, factor_txt_dict_save_path)
            save(total, factor_txt_total_save_path)

save(dicts, factor_txt_dict_save_path)
save(total, factor_txt_total_save_path)

print(total)
factors = list(dicts[-1].keys())
counts = list(dicts[-1].values())

sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i])

sorted_factors = [factors[i] for i in sorted_indices]
sorted_counts = [counts[i] for i in sorted_indices]

with open(factor_csv_save_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([''] + [f'{classification_dict[i]}' for i in range(len(dicts))])
    for key in sorted_factors:
        row = [key] + [d.get(key, 0) / total[i] for i, d in enumerate(dicts)]
        writer.writerow(row)

# save_path_ids = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/data_analysis_jaspar/ids.txt"
# save_path_matrixes = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/data_analysis_jaspar/matrixes.txt"
# url = "http://jaspar.genereg.net/api/v1/matrix/?tax_group=vertebrates"
# validated_matrix_ids = []
# matrixes = []

# def save(responses, factor_save_path):
#     with open(factor_save_path, 'w') as f:
#         for item in responses:
#             f.write("%s/n" % item)

# while url:
#     response = requests.get(url)
#     data = json.loads(response.text)
#     validated_matrices = [matrix for matrix in data['results'] if matrix['collection'] == 'CORE']
#     ids = [matrix['matrix_id'] for matrix in validated_matrices]
#     validated_matrix_ids.extend(ids)
#     url = data['next']
#     print(url)

# save(validated_matrix_ids, save_path_ids)

# for matrix_id in validated_matrix_ids:
#     matrix_url = f"http://jaspar.genereg.net/api/v1/matrix/{matrix_id}/"
#     matrix_response = requests.get(matrix_url)
#     matrix_data = json.loads(matrix_response.text)
#     matrixes.append(matrix_data)

# save(matrixes, save_path_matrixes)