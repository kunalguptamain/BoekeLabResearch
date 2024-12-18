
import sys
sys.path.insert(0, r'C:/Users/kunal/Documents/BoekeLabResearch/')
import methods.general as general
from pyjaspar import jaspardb
import time

promoter_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_old_model/compare_sequences_8.txt"
factor_dict_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_old_model/factors_old_dict_8.txt"
factor_total_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_old_model/factors_old_total_8.txt"
factor_csv_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_old_model/factors_old_8.csv"

label_dict = {
    0: "HUMAN_CODING",
    1: "HUMAN_NON_CODING",
    2: "MOUSE_CODING",
    3: "MOUSE_NON_CODING",
    4: "NON_PROMOTER",
    5: "HUMAN_CODING_GENEREATED",
}

threshold_dissimilarity = 3

jdb_obj = jaspardb(release='JASPAR2024')

motifs = jdb_obj.fetch_motifs(
    collection = ['CORE'],
    all_versions = False,
    tax_group='Vertebrates',
)
start_time = time.time()
if __name__ == "__main__":
    general.analyse_intake_data(threshold_dissimilarity=threshold_dissimilarity,
                                motifs=motifs,
                                label_dict=label_dict,
                                promoter_path=promoter_path,
                                factor_dict_save_path=factor_dict_save_path,
                                factor_total_save_path=factor_total_save_path,
                                factor_csv_save_path=factor_csv_save_path)
end_time = time.time()
print((start_time - end_time))