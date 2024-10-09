from standardized_dna_diffusion_pipeline.sequence_set import SequenceSet
from standardized_dna_diffusion_pipeline.analyzer import Analyzer
from pyjaspar import jaspardb
from functools import partial
from Bio import motifs
from Bio.Seq import Seq

data_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/pre_generation_sequence_set.json'
sequences = SequenceSet.from_file(data_path)
print("loaded data")

analyzer = Analyzer(sequences)
coding_hk_query = analyzer.query_by_class("CODING_HK_DNA")
noncoding_nonhk_query = analyzer.query_by_class("NON_CODING_NON_HK_DNA")
junk_dna_query = analyzer.query_by_class("JUNK_DNA")
dual_query = coding_hk_query + junk_dna_query + noncoding_nonhk_query

analyzer.perform_analysis(
    analyzer.analyze_gc_content,
    dual_query
)

jaspar_db = jaspardb(release="JASPAR2024")
tbp_factors = jaspar_db.fetch_motif_by_id("MA0108.1")
analyzer.perform_analysis(
    analysis_func=analyzer.analysis_tfs,
    sequences=dual_query,
    params={
        "motifs": [tbp_factors],
        "threshold": 0.8,
        "name": "TBP Binding Site"
    }
)

print("finished analyzing")

analyzer.plot_analysis(
    [junk_dna_query, coding_hk_query, noncoding_nonhk_query],
    ["TBP Binding Site"],
    ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA"],
    normalize=False
)
