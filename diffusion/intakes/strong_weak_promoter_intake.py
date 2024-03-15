from Bio.Seq import Seq
import numpy as np
import os

one_hot_encode_mapping = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1],
}

classification_dict = {
    "STRONG": 2,
    "WEAK": 1,
    "NON": 0,
}

promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/strong_weak_promoters/promoter.txt"
non_promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/strong_weak_promoters/non_promoter.txt"
numpy_abs_path = "C:/Users/kunal/Documents/npy_data/"
     
def one_hot_encode(seq, total_len):
    encoded = [one_hot_encode_mapping[i] for i in seq]
    length = len(encoded)
    for _ in range(total_len - length): encoded.append([0, 0, 0, 0])
    return np.array(encoded)

save_data = []
save_labels = []

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
            if(not placed): continue
            sequence = Seq(line.strip())
            if(not forward):sequence = sequence.reverse_complement()
            save_data.append(one_hot_encode(str(sequence), 88).T)
            save_labels.append(classification)

print(save_data[0][:, :10])

with open(non_promoter_abs_path, 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        if(i == 1632): break
        save_data.append(one_hot_encode(line.strip(), 88).T)
        save_labels.append(0)

save_data = np.array(save_data)
save_labels = np.array(save_labels)
np.save(os.path.join(numpy_abs_path, "save_data_promoter_strong_weak.npy"), save_data)
np.save(os.path.join(numpy_abs_path, "save_labels_promoter_strong_weak.npy"), save_labels)

print(save_data.shape)
print(save_labels.shape)