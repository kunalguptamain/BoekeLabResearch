import sys
sys.path.insert(0, r'C:/Users/kunal/Documents/BoekeLabResearch/')

from pyjaspar import jaspardb
import numpy as np
import csv
from Bio.Seq import Seq
import json
import methods.general as general

promoter_abs_path = r"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/enhancer_generation/data/compare_sequences.txt"
factor_txt_dict_save_path = r"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/enhancer_generation/train_data_analysis/factors_enhancers_dict.txt"
factor_txt_total_save_path = r"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/enhancer_generation/train_data_analysis/factors_enhancers_total.txt"
factor_csv_save_path = r"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/enhancer_generation/train_data_analysis/factors_enhancers.csv"

label_dict = {
    0: "0 -- 1.93483998",
    1: "3.22714918",
    2: "4.33691651",
    3: "6.96985935",
    4: "70",
    5: "GENERAL"
}

# jdb_obj = jaspardb(release='JASPAR2024')

# motifs = jdb_obj.fetch_motifs(
#     collection = ['CORE'],
#     all_versions = False,
#     tax_group='Vertebrates',
# )

# dicts = [{}, {}, {}, {}, {}, {}]
# total = [0, 0, 0, 0, 0, 0]

# with open(promoter_abs_path) as file:
#     sequence_set = json.load(file)

# cycle = 100
# i = 0
# for sequence, label in sequence_set:
#     print(i)
#     general.get_and_sort_motifs(sequence, 2, motifs, [dicts[label], dicts[-1]])
#     total[label] += 1
#     total[-1] += 1
#     i += 1
#     if(i % cycle == 0):
#         general.save_obj(dicts, factor_txt_dict_save_path)
#         general.save_obj(total, factor_txt_total_save_path)

# general.save_obj(dicts, factor_txt_dict_save_path)
# general.save_obj(total, factor_txt_total_save_path)

with open(factor_txt_dict_save_path) as file:
    dicts = json.load(file)

with open(factor_txt_total_save_path) as file:
    total = json.load(file)

general.generate_motif_analysis_csv(dicts, total, label_dict, factor_csv_save_path)
