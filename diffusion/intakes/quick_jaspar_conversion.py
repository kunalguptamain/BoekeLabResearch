import numpy as np
import os

num_files = 1704
height = 4
width = 40
trials_per_site = 20
binding_sites_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/chipseq_binding_sites"

def read_jaspar(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        matrix = [list(map(int, line[line.find('[')+1:line.find(']') -1].split())) for line in lines]
    return np.array(matrix, dtype=int)

def normalize(matrix):
    # Normalize the matrix so that the sum of each column is 1
    column_sums = matrix.sum(axis=0)
    return matrix[:, :] / column_sums[np.newaxis, :]

# Replace 'your_directory' with the path to your directory of .jaspar files
for i, file_name in enumerate(os.listdir(binding_sites_abs_path)):
    if file_name.endswith('.jaspar'):
        file_path = os.path.join(binding_sites_abs_path, file_name)
        matrix = read_jaspar(file_path)
        normalized_matrix = normalize(matrix)
        print(file_name, normalized_matrix.tolist())