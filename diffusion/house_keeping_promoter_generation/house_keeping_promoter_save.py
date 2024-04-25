import numpy as np
from enum import Enum
import os
import json

class Class_Labels(Enum):
    HUMAN_CODING = 0
    HUMAN_NON_CODING = 1
    MOUSE_CODING = 2
    MOUSE_NON_CODING = 3
    NON_PROMOTER = 4

file_name_to_class_dict = {
    "hg38_EYl13.fa.txt": Class_Labels.HUMAN_CODING.value, #human
    "hg38_YHek6.fa.txt": Class_Labels.HUMAN_NON_CODING.value, #human non coding
    "mm10_8iiKK.fa.txt": Class_Labels.MOUSE_CODING.value, #house mouse
    "mm10_R33Vk.fa.txt": Class_Labels.MOUSE_NON_CODING.value, #house mouse non coding
}

promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data/"
non_promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data/nonpromoter_upper.txt"
select_promoters_human_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data/Housekeeping_TranscriptsHuman.csv"
select_promoters_mouse_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data/Housekeeping_TranscriptsMouse.csv"
save_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data/compare_sequences_8.txt"

def read_sequences(file_path):
     with open(file_path, 'r') as file:
        lines = file.readlines()
        sequences = []
        sequence = ''
        visited = False
        for i in range(len(lines)):
            if i % 11 == 0:
                if(visited):
                    sequences.append((sequence_name, sequence))
                    sequence = ''
                visited = True
                tokens = lines[i].split()
                sequence_name = tokens[1]
                continue
            sequence += lines[i].rstrip()
        return sequences
    
gene_name_to_sequences = {}

for i, file_name in enumerate(os.listdir(promoter_abs_path)):
    if not file_name.endswith('.txt'): continue
    if file_name not in file_name_to_class_dict: continue
    file_label = file_name_to_class_dict[file_name]
    file_path = os.path.join(promoter_abs_path, file_name)
    sequences = read_sequences(file_path)
    for sequence_data in sequences:
        gene_name, sequence = sequence_data
        simplified_gene_name = gene_name[0:gene_name.index('_')]
        if simplified_gene_name not in gene_name_to_sequences:
            gene_name_to_sequences[simplified_gene_name] = []
        gene_name_to_sequences[simplified_gene_name].append((sequence, file_label))

save_data = []

def readCSV(file_path):
    with open(file_path) as file:
        lines = file.readlines()
        for line in lines:
            gene_symbol = line.split(';')[1]
            if gene_symbol not in gene_name_to_sequences: continue
            relevant_sequences = gene_name_to_sequences[gene_symbol]
            for sequence_data in relevant_sequences:
                sequence, label = sequence_data
                sequence = sequence[300:-50]
                save_data.append((sequence, label))

readCSV(select_promoters_human_path)
readCSV(select_promoters_mouse_path)

for sequence_set in gene_name_to_sequences.values():
    for sequence, label in sequence_set:
        if(label not in [Class_Labels.HUMAN_NON_CODING.value, \
                         Class_Labels.MOUSE_NON_CODING.value]): continue
        sequence = sequence[300:-50]
        save_data.append((sequence, label))

data_collected = len(save_data)

with open(non_promoter_abs_path) as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        if(i == data_collected): break
        sequence = line.strip()[25:-25]
        save_data.append((sequence, Class_Labels.NON_PROMOTER.value))

with open(save_path, 'w') as file:
    file.write(json.dumps(save_data))
