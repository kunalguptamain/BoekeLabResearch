import sys
sys.path.insert(0, r'C:/Users/kunal/Documents/BoekeLabResearch/')

import methods.general as general
import pandas as pd

ENHANCER_LOAD_PATH = r'C:/Users/kunal/Documents/BoekeLabResearch/diffusion/enhancer_generation/data/Peng Gen Res 2020.csv'
SAVE_PATH = r"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/enhancer_generation/data/enhancers.txt"

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
for sequence in sequences:
    if("N" in sequence): continue
    label = get_bucket(ENHANCER_SCORE_BUCKETS, enhancer_scores[i])
    save_data.append((sequence, label))

general.save_obj(save_data, SAVE_PATH)