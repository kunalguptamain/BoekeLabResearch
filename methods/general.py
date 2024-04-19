from pyjaspar import jaspardb
import numpy as np
import csv
from Bio.Seq import Seq
import json
import copy
import pandas as pd

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
        totals,
        class_label_dict,
        csv_save_path,
        sort_label = None,
        existing_csv_path = None,
    ):
    #creating consolidated columns
    net_labels = list(class_label_dict.values())
    if existing_csv_path:
        df = pd.read_csv(existing_csv_path)
        net_labels = net_labels + list(df.columns)
    net_labels = list(set([label for label in net_labels if label not in ["GENERAL", "Unnamed: 0"]]))
    net_labels.sort()
    net_labels = ["FACTOR"] + net_labels

    #creating consolidated rows
    total_factors = set()
    for d in dicts:
        for key in d.keys():
            total_factors.add(key)
    if existing_csv_path:
        for factor in list(df.iloc[:, 0]):total_factors.add(factor)
    total_factors = list(total_factors)

    #sorting factors
    if sort_label in class_label_dict.values():
        dict_index = list(class_label_dict.values()).index(sort_label)
        total_factors.sort(key=lambda x: dicts[dict_index].get(x, 0), reverse=True)
    elif existing_csv_path and sort_label in list(df.columns):
        total_factors.sort(key=lambda x: df[df.iloc[:, 0] == x][sort_label].values[0] if x in df.iloc[:, 0].values else 0, reverse=True)
    
    #adding values in 
    if existing_csv_path:
        df.set_index(df.columns[0], inplace=True)
    save_frame = pd.DataFrame(columns=net_labels)
    save_frame["FACTOR"] = total_factors
    for i, col_name in enumerate(net_labels):
        if(col_name == "FACTOR"): continue
        for j, factor in enumerate(total_factors):
            if col_name in class_label_dict.values():
                dict_index = list(class_label_dict.values()).index(col_name)
                save_frame.loc[j, col_name] = dicts[dict_index][factor] / totals[dict_index] if factor in dicts[dict_index] else 0
            elif existing_csv_path:
                save_frame.loc[j, col_name] = df.loc[factor, col_name] if (factor in df.index) and (col_name in df.columns) else 0
    
    save_frame.to_csv(csv_save_path, index=False)
        

def json_load(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)
    
def analyse_intake_data(threshold_dissimilarity,
    motifs,
    label_dict,
    promoter_path,
    factor_dict_save_path = "factor_dict.txt",
    factor_total_save_path = "factor_total.txt",
    factor_csv_save_path = "factor_csv.csv",
    sort_label = None,
    existing_csv_path = None,
):
    dicts = [{} for _ in range(len(label_dict))]
    total = [0 for _ in range(len(label_dict))]

    sequence_set = json_load(promoter_path)

    i = 0
    for sequence, label in sequence_set:
        i += 1
        if(i == 20): break
        print(i)
        if(label not in label_dict): continue
        get_and_sort_motifs(sequence, threshold_dissimilarity, motifs, [dicts[label]])
        total[label] += 1

    save_obj(dicts, factor_dict_save_path)
    save_obj(total, factor_total_save_path)

    generate_motif_analysis_csv(dicts, total, label_dict, factor_csv_save_path, sort_label, existing_csv_path)