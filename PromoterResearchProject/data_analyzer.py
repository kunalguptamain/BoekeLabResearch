from standardized_dna_diffusion_pipeline.sequence_set import SequenceSet
from standardized_dna_diffusion_pipeline.analyzer import Analyzer, FrequentMotifs
from pyjaspar import jaspardb
from functools import partial
from typing import List

data_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/sequences_set_states/32-128post_generation_sequence_set_3_and_2.json'
sequences = SequenceSet.from_file(data_path)
print("loaded data")

analyzer = Analyzer(sequences)

generated_coding_hk_query = analyzer.query_by_class_and_generated("CODING_HK_DNA", True)
coding_hk_query = analyzer.query_by_class_and_generated("CODING_HK_DNA", False)
noncoding_nonhk_query = analyzer.query_by_class("NON_CODING_NON_HK_DNA")
generated_coding_nonhk_query = analyzer.query_by_class_and_generated("CODING_NON_HK_DNA", True)
coding_nonhk_query = analyzer.query_by_class_and_generated("CODING_NON_HK_DNA", False)
junk_dna_query = analyzer.query_by_class("JUNK_DNA")

dual_query = coding_hk_query + junk_dna_query + noncoding_nonhk_query + coding_nonhk_query + generated_coding_hk_query + generated_coding_nonhk_query

analyzer.perform_analysis(
    analyzer.analyze_gc_content,
    dual_query
)

jaspar_db = jaspardb(release="JASPAR2024")
gabpa_factor = jaspar_db.fetch_motif_by_id("MA0062.2")
analyzer.perform_analysis(
    analysis_func=partial(
        analyzer.analysis_tfs,
        motifs = [gabpa_factor],
        threshold = 0.8,
        input_name = "GABPA"
    ),
    sequences=dual_query,
)

# jaspar_db = jaspardb(release="JASPAR2024")
# list = ["MA0158.1", "MA0594.1", "MA0650.3", "MA0650.4", "MA0911.1", "MA0911.2"]
# hoxa_sites = [jaspar_db.fetch_motif_by_id(string) for string in list]
# analyzer.perform_analysis(
#     analysis_func=partial(
#         analyzer.analysis_tfs,
#         motifs = hoxa_sites,
#         threshold = 0.8,
#         input_name = "HOX"
#     ),
#     sequences=dual_query,
# )
# list = ["MA0151.1", "MA0601.1", "MA0601.2", "MA0602.1", "MA0602.2",]
# arid_sites = [jaspar_db.fetch_motif_by_id(string) for string in list]
# analyzer.perform_analysis(
#     analysis_func=partial(
#         analyzer.analysis_tfs,
#         motifs = arid_sites,
#         threshold = 0.8,
#         input_name = "ARID"
#     ),
#     sequences=dual_query,
# )
# tbp_factor = jaspar_db.fetch_motif_by_id("MA0108.1")
# analyzer.perform_analysis(
#     analysis_func=partial(
#         analyzer.analysis_tfs,
#         motifs = [tbp_factor],
#         threshold = 0.8,
#         input_name = "TBP Binding Site"
#     ),
#     sequences=dual_query,
# )

# all_factors = jaspar_db.fetch_motifs(
#     collection = ['CORE'],
#     all_versions = False,
#     tax_group='Vertebrates',
#     species=10090
# )

# analyzer.perform_analysis(
#     analysis_func=partial(
#         analyzer.analysis_tfs,
#         motifs = all_factors,
#         threshold = 0.8,
#         input_name = "All Vertebrate Factors"
#     ),
#     sequences=dual_query,
# )

#-------------------------------------------------------------------
motif_set = FrequentMotifs()
motifs = {
    "TATA": FrequentMotifs.tata_motif,
    "Initiator": FrequentMotifs.initiator_motif,
    "CCAAT": FrequentMotifs.ccaat_motif,
    "GC BOX": FrequentMotifs.gc_motif,
}

for name, motif in motifs.items():
    analyzer.perform_analysis(
        analysis_func=partial(
            analyzer.analysis_tfs,
            motifs = [motif],
            threshold = 0.8,
            input_name = name
        ),
        sequences=dual_query,
    )

# # pou5f1_sox2 = jaspar_db.fetch_motif_by_id("MA0142.1")
# analyzer.perform_analysis(
#     analysis_func=partial(
#         analyzer.analysis_tfs,
#         motifs = [pou5f1_sox2],
#         threshold = 0.8,
#         input_name = "pou5f1::sox2"
#     ),
#     sequences=dual_query,
# )
# housekeeping_factors = [
#     jaspar_db.fetch_motif_by_id("MA1535.1")
# ]
# housekeeping_factors = [factor for factor in housekeeping_factors if factor]
# analyzer.perform_analysis(
#     analysis_func=partial(
#         analyzer.analysis_tfs,
#         motifs = housekeeping_factors,
#         threshold = 0.8,
#         input_name = "NR2C1"
#     ),
#     sequences=dual_query,
# )

# fm = FrequentMotifs()
# regulartory_list = {
#     "DRE": fm.dre,
#     "BREU": fm.breu,
#     "RBM": fm.rbm
# }
# print(fm.dre.pssm)
# for name, motif in regulartory_list.items():
#     analyzer.perform_analysis(
#         analysis_func=partial(
#             analyzer.analysis_tfs,
#             motifs = [motif],
#             threshold = 0.8,
#             input_name = name,
#         ),
#         sequences=dual_query,
#     )


# print("finished analyzing")

# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, noncoding_nonhk_query, coding_nonhk_query],
#     ["GC Content"],
#     ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA", "CODING NON HK DNA"],
#     normalize=False,
#     reduce_noise=True
# )

# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, noncoding_nonhk_query, coding_nonhk_query],
#     ["GABPA"],
#     ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA", "CODING NON HK DNA"],
#     normalize=False,
#     reduce_noise=True
# )

# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, noncoding_nonhk_query, coding_nonhk_query],
#     ["Initiator"],
#     ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA", "CODING NON HK DNA"],
#     normalize=False,
#     reduce_noise=True
# )

# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, noncoding_nonhk_query, coding_nonhk_query],
#     ["TATA"],
#     ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA", "CODING NON HK DNA"],
#     normalize=False,
#     reduce_noise=True
# )

# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, noncoding_nonhk_query, coding_nonhk_query],
#     ["pou5f1::sox2", "TBP Binding Site"],
#     ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA", "CODING NON HK DNA"],
#     normalize=False,
#     reduce_noise=True,
# )


# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, generated_coding_hk_query],
#     ["GABPA"],
#     ["JUNK DNA", "CODING HK DNA", "GENERATED CODING HK DNA"],
#     normalize=False,
#     reduce_noise=False,
# )


# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, generated_coding_hk_query],
#     ["GC Content"],
#     ["JUNK DNA", "CODING HK DNA", "GENERATED CODING HK DNA"],
#     normalize=False,
#     reduce_noise=True,
# )


# for name in motifs.keys():
#     analyzer.plot_analysis(
#         [junk_dna_query, coding_hk_query, generated_coding_hk_query],
#         [name],
#         ["JUNK DNA", "CODING HK DNA", "GENERATED CODING HK DNA"],
#         normalize=False,
#         reduce_noise=True
#     )

plotValuesList = {
    "nonhk": [junk_dna_query, coding_nonhk_query, generated_coding_nonhk_query],
    "hk": [junk_dna_query, coding_hk_query, generated_coding_hk_query],
}
labelsList = {
    "nonhk": ["JUNK DNA", "CODING NON HK DNA", "GENERATED CODING NON HK DNA"],
    "hk": ["JUNK DNA", "CODING HK DNA", "GENERATED CODING HK DNA"]
}

plotValues = plotValuesList["hk"]
labels = labelsList["hk"]

analyzer.plot_analysis(
    plotValues,
    ["GABPA"],
    labels,
    normalize=False,
    reduce_noise=True,
)
analyzer.plot_analysis(
    plotValues,
    ["GC Content"],
    labels,
    normalize=False,
    reduce_noise=True,
)
for name in motifs.keys():
    analyzer.plot_analysis(
        plotValues,
        [name],
        labels,
        normalize=False,
        reduce_noise=True
    )


# analyzer.plot_analysis(
#     plotValues,
#     ["GC Content"],
#     labels,
#     normalize=False,
#     reduce_noise=True,
# )


# for name in motifs.keys():
#     analyzer.plot_analysis(
#         plotValues,
#         [name],
#         labels,
#         normalize=False,
#         reduce_noise=True
#     )



# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, noncoding_nonhk_query, coding_nonhk_query],
#     ["All Vertebrate Factors"],
#     ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA", "CODING NON HK DNA"],
#     normalize=False
# )

# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, noncoding_nonhk_query, coding_nonhk_query],
#     ["TBP Binding Site"],
#     ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA", "CODING NON HK DNA"],
#     normalize=False
# )

# analyzer.plot_analysis(
#     [junk_dna_query, coding_hk_query, noncoding_nonhk_query, coding_nonhk_query],
#     ["GC Content"],
#     ["JUNK DNA", "CODING HK DNA", "NON CODING NON HK DNA", "CODING NON HK DNA"],
#     normalize=False
# )
