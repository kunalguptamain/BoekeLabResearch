from pyjaspar import jaspardb
import numpy as np
import csv
from Bio.Seq import Seq
import json
import copy

def one_hot_encode(seq, total_len):
    one_hot_encode_mapping = {
        'A': [1, 0, 0, 0],
        'C': [0, 1, 0, 0],
        'G': [0, 0, 1, 0],
        'T': [0, 0, 0, 1],
    }
    encoded = [one_hot_encode_mapping[i] for i in seq]
    length = len(encoded)
    for _ in range(total_len - length): encoded.append([0, 0, 0, 0])
    return np.array(encoded).T


def count_single_motif(sequence, target_motif, threshold):
    pssm = target_motif.pwm.log_odds()
    search = list(pssm.search(Seq(sequence), threshold=threshold))
    return len(search)

def get_motifs(sequence, threshold_dissimilarity, motifs):
    '''
    Given a sequence, the method will find the location of motifs within the threshold
    disimilarity of the Posisition Weight Matrix. It will count them up and add them to a 
    dictionary where key is MOTIF NAME and value is 
    '''
    return_dict = {}
    for i, motif in enumerate(motifs):
        count = count_single_motif(sequence, motif, threshold_dissimilarity)
        if(count == 0): continue
        if motif.name not in return_dict: return_dict[motif.name] = 0
        return_dict[motif.name] += count
    return return_dict
        
def get_and_sort_motifs(sequence, threshold_dissimilarity, motifs, sort_into_dicts):
    '''
    Adds motifs collected from getmotifs() into a list of dictionaries
    '''
    for motif, count in get_motifs(sequence, threshold_dissimilarity, motifs).items():
        for dict_r in sort_into_dicts:
            if motif not in dict_r: dict_r[motif] = 0
            dict_r[motif] += count

def save_obj(object, path):
    with open(path, 'w') as file:
        file.write(json.dumps(object))

def generate_motif_analysis_csv(
        dicts,
        total,
        class_label_dict,
        csv_save_path,
    ):

    factors = list(dicts[-1].keys())
    counts = list(dicts[-1].values())

    sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i])
    sorted_factors = [factors[i] for i in sorted_indices]

    with open(csv_save_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([''] + [f'{class_label_dict[i]}' for i in range(len(dicts))])
        for key in sorted_factors:
            row = [key] + [d.get(key, 0) / total[i] for i, d in enumerate(dicts)]
            writer.writerow(row)
