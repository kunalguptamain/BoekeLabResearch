from pyjaspar import jaspardb
import numpy as np
import csv
from Bio.Seq import Seq
from multiprocessing import Pool, Manager, current_process
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
    return np.array(encoded).T.tolist()