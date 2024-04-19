
import sys
sys.path.insert(0, r'C:/Users/kunal/Documents/BoekeLabResearch/')
import methods.general as general
from pyjaspar import jaspardb

promoter_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_old_model/compare_sequences.txt"
factor_dict_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_old_model/factors_old_dict.txt"
factor_total_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_old_model/factors_old_total.txt"
factor_csv_save_path = f"C:/Users/kunal/Documents/BoekeLabResearch/diffusion/house_keeping_promoter_generation/data_analysis_old_model/factors_old.csv"

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

general.analyse_intake_data(threshold_dissimilarity=threshold_dissimilarity,
                            motifs=motifs,
                            label_dict=label_dict,
                            promoter_path=promoter_path,
                            factor_dict_save_path=factor_dict_save_path,
                            factor_total_save_path=factor_total_save_path,
                            factor_csv_save_path=factor_csv_save_path)