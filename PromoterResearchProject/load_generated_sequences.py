from standardized_dna_diffusion_pipeline.sequence_set import SequenceSet

data_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/sequences_set_states/pre_generation_sequence_set_limited.json'
new_data_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/sequences_set_states/128_1x2x2x2x2x4x4x8x8x16_post_generation_sequence_set_3_and_2.json'

generated_promoters_path_1 = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/generated_sequences/generated_class_3_promoters_128_1x2x2x2x2x4x4x8x8x16.txt'
generated_promoters_path_2 = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/generated_sequences/generated_class_2_promoters_128_1x2x2x2x2x4x4x8x8x16.txt'
sequences = SequenceSet.from_file(data_path)

with open(generated_promoters_path_1, 'r') as file:
    lines = file.readlines()
    for line in lines:
        sequences.add_sequence(
            "CODING_HK_DNA",
            line.strip(),
            is_generated=True
        )


with open(generated_promoters_path_2, 'r') as file:
    lines = file.readlines()
    for line in lines:
        sequences.add_sequence(
            "CODING_NON_HK_DNA",
            line.strip(),
            is_generated=True
        )

sequences.save(new_data_path)
