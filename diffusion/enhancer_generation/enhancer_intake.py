import sys
sys.path.insert(0, r'C:/Users/kunal/Documents/BoekeLabResearch/')

import methods.general as general
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import os

ENHANCER_LOAD_PATH = r'C:/Users/kunal/Documents/BoekeLabResearch/diffusion/enhancer_generation/data/Peng Gen Res 2020.csv'
NUMPY_ABS_PATH = r"C:/Users/kunal/Documents/npy_data/"

df = pd.read_csv(ENHANCER_LOAD_PATH)
sequence_numbers = df['Sequence Number'].to_numpy()
sequences = df['Sequence'].to_numpy()
chr_positions = df['chr_pos'].to_numpy()
enhancer_scores = df['Enhancer Score'].to_numpy()

def get_bucket(thresholds, enhancer_score):
    for i, threshold in enumerate(thresholds):
        if(enhancer_score <= threshold):
            return i
        
ENHANCER_SCORE_BUCKETS = [1.93483998,  3.22714918,  4.33691651,  6.96985935, 70]
MAX_SEQ_LEN = 9200

save_data = []
save_labels = []
for i, sequence in enumerate(sequences):
    print(i)
    if("N" in sequence): continue
    reversed_sequence = str(Seq(sequence).reverse_complement())
    label = get_bucket(ENHANCER_SCORE_BUCKETS, enhancer_scores[i])
    save_data.append(general.one_hot_encode(sequence, MAX_SEQ_LEN))
    save_data.append(general.one_hot_encode(reversed_sequence, MAX_SEQ_LEN))
    save_labels.append(label)
    save_labels.append(label)


save_data = np.array(save_data)
save_labels = np.array(save_labels)
print(save_data.shape)
print(save_labels.shape)
np.save(os.path.join(NUMPY_ABS_PATH, "save_data_enhancer.npy"), save_data)
np.save(os.path.join(NUMPY_ABS_PATH, "save_labels_enhancer.npy"), save_labels)