import numpy as np
import os

file_name_to_class_dict = {
    "hg38_EYl13.fa.txt": 1, #human
    "hg38_YHek6.fa.txt": 1, #human non coding
    "mm10_8iiKK.fa.txt": 2, #house mouse
    "mm10_R33Vk.fa.txt": 2, #house mouse non coding
    "rheMac8_BhvOq.fa.txt": 3, #monke ooh ah bannana
    "galGal5_liMv6.fa.txt": 4, #birb :>
    "rn6_idvdz.fa.txt": 5, #rat
    "canFam3_EhUhY.fa.txt": 6#dog
}

one_hot_encode_mapping = {
    'A': [1, 0, 0, 0, 0, 0, 0, 0],
    'C': [0, 1, 0, 0, 0, 0, 0, 0],
    'G': [0, 0, 1, 0, 0, 0, 0, 0],
    'T': [0, 0, 0, 1, 0, 0, 0, 0],
    'N': [0, 0, 0, 0, 1, 0, 0, 0]
}

promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/promoters"
numpy_abs_path = "C:/Users/kunal/Documents/npy_data/"

def read_sequences(file_path):
     with open(file_path, 'r') as file:
        lines = file.readlines()
        sequences = []
        sequence = ''
        for i in range(len(lines)):
            if i % 11 == 0:
                sequences.append(sequence)
                sequence = ''
                continue
            sequence += lines[i].rstrip()
        return sequences
     
def one_hot_encode(seq):
    return np.array([one_hot_encode_mapping[i] for i in seq])

save_data = []
save_labels = []

for i, file_name in enumerate(os.listdir(promoter_abs_path)):
    if not file_name.endswith('.txt'): continue
    label = file_name_to_class_dict[file_name]
    file_path = os.path.join(promoter_abs_path, file_name)
    sequences = read_sequences(file_path)
    for sequence in sequences:
        encoding = one_hot_encode(sequence)
        if encoding.shape != (600, 4):continue
        save_data.append(encoding)
        save_labels.append(label)

save_data = np.array(save_data)
save_labels = np.array(save_labels)
np.save(os.path.join(numpy_abs_path, "save_data_promoter_alt_mapping.npy"), save_data)
np.save(os.path.join(numpy_abs_path, "save_labels_promoter_alt_mapping.npy"), save_labels)
