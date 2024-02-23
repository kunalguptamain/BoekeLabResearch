import numpy as np
import os

num_files = 1704
height = 4
width = 40
trials_per_site = 20
binding_sites_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/bindingsites"


def read_jaspar(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        matrix = [list(map(int, line[line.find('[')+1:line.find(']') -1].split())) for line in lines]
    return np.array(matrix, dtype=int)

def normalize(matrix):
    # Normalize the matrix so that the sum of each column is 1
    column_sums = matrix.sum(axis=0)
    return matrix[:, :] / column_sums[np.newaxis, :]

def generate_multinomial_outcomes(n, matrix):
    # Apply the multinomial distribution to each column
    return np.apply_along_axis(lambda col: np.random.multinomial(n, col), axis=0, arr=matrix)

def padding(matrix):
    return np.pad(matrix, ((0, 0), (0, width - matrix.shape[1])), constant_values=0)

def generate_sequences(matrix):
    sequence_matrix = np.zeros((trials_per_site, width))
    for col in range(matrix.shape[1]):
        sum_row = 0
        for row in range(matrix.shape[0]):
            sequence_matrix[sum_row : sum_row + matrix[row][col], col] = row + 1
            sum_row += matrix[row][col]

    np.apply_along_axis(lambda x: np.random.shuffle(x), 0, sequence_matrix)
    sequence_matrix = sequence_matrix.astype(np.int16)

    one_hot_matrix = np.zeros((trials_per_site, height, width))
    for i in range(trials_per_site):
        one_hot_encoding = np.squeeze(np.eye(5)[sequence_matrix[i, :].reshape(-1)])
        one_hot_encoding = np.delete(one_hot_encoding, [0], axis=1)
        one_hot_encoding = np.matrix.transpose(one_hot_encoding)
        one_hot_encoding = np.expand_dims(one_hot_encoding, axis=0)
        one_hot_matrix[i, :, :] = np.expand_dims(one_hot_encoding, axis=0)

    return one_hot_matrix

save_data = np.zeros((num_files * trials_per_site, height, width))
save_labels = np.zeros((num_files * trials_per_site))

# Replace 'your_directory' with the path to your directory of .jaspar files
for i, file_name in enumerate(os.listdir(binding_sites_abs_path)):
    if file_name.endswith('.jaspar'):
        file_path = os.path.join(binding_sites_abs_path, file_name)
        matrix = read_jaspar(file_path)
        normalized_matrix = normalize(matrix)
        outcomes = generate_multinomial_outcomes(trials_per_site, normalized_matrix)
        padded_outcomes = padding(outcomes)
        sequences = generate_sequences(padded_outcomes)
        trial_sums = trials_per_site * i
        save_data[trial_sums:(trial_sums + trials_per_site), :, :] = sequences
        save_labels[trial_sums:(trial_sums + trials_per_site)] = i

numpy_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/"
np.save(os.path.join(numpy_abs_path, "save_data.npy"), save_data)
np.save(os.path.join(numpy_abs_path, "save_labels.npy"), save_labels)
