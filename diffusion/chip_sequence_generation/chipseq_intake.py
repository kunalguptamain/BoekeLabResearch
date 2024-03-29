import numpy as np
import os

file_name_to_class_dict = {
    "GATA1.fa": 1,
    "JUND.fa": 2,
    "MYC.fa": 3,
    "POLR2A.fa": 4,
    "TBP.fa": 5,
    "TCF12.fa": 6,
    "USF1.fa": 7
}
letters_to_numbers = {
    "A": 1,
    "C": 2,
    "G": 3,
    "T": 4,
    "N": 5
}
chipseq_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/chipseq/"
numpy_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/"
max_sequence_length = 1000
height = 4

def read_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        file_sequences = []
        for i, line in enumerate(lines):
            if(i % 2 == 0): continue
            cleaned_sequence = line[:-2].upper()
            if(len(cleaned_sequence) > max_sequence_length): continue
            # Create a numpy array of the sequence
            sequence_array = np.array(list(cleaned_sequence))
            # Convert the sequence using dictionary mapping
            converted_sequence = np.array([letters_to_numbers[letter] for letter in sequence_array])

            one_hot_encoding = np.squeeze(np.eye(5)[converted_sequence]).astype('int8')
            one_hot_encoding = np.delete(one_hot_encoding, [0], axis=1)
            one_hot_encoding = np.matrix.transpose(one_hot_encoding)
            padding = ((0, 0), (0, max(0, max_sequence_length - one_hot_encoding.shape[1])))
            one_hot_encoding = np.pad(one_hot_encoding, padding)
            if one_hot_encoding.shape[1] != max_sequence_length: print(one_hot_encoding.shape)

            file_sequences.append(one_hot_encoding)
    return file_sequences

save_data = []
save_labels = []

for i, file_name in enumerate(os.listdir(chipseq_abs_path)):
    if file_name.endswith('.fa'):
        file_path = os.path.join(chipseq_abs_path, file_name)
        gen_sequences = read_file(file_path)
        for i in range(len(gen_sequences)): 
            save_data.append(gen_sequences[i])
            save_labels.append(file_name_to_class_dict[file_name])

save_data = np.array(save_data)
save_labels = np.array(save_labels)
np.save(os.path.join(numpy_abs_path, "save_data_chipseq.npy"), save_data)
np.save(os.path.join(numpy_abs_path, "save_labels_chipseq.npy"), save_labels)
