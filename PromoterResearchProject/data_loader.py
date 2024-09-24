from standardized_dna_diffusion_pipeline.sequence_set import SequenceSet
from enum import Enum
from Bio import SeqIO

mouse_junk_dna_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/mouseDNAData/junk_dna.txt'
mouse_coding_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/mouseDNAData/mouse_coding.txt'
mouse_non_coding_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/mouseDNAData/mouse_non_coding.txt'
house_keeping_list_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/mouseDNAData/house_keeping.txt'

class PromoterTypes(str, Enum):
    JUNK_DNA = "JUNK_DNA"
    NON_CODING_NON_HK_DNA = "NON_CODING_NON_HK_DNA"
    CODING_NON_HK_DNA = "CODING_NON_HK_DNA"
    CODING_HK_DNA = "CODING_HK_DNA"

sequence_set = SequenceSet(
    dna_class_types=[
        PromoterTypes.JUNK_DNA.value,
        PromoterTypes.NON_CODING_NON_HK_DNA.value,
        PromoterTypes.CODING_NON_HK_DNA.value,
        PromoterTypes.CODING_HK_DNA.value,
    ],
    sequence_padded_length=1024
)

house_keeping_list = []
with open(house_keeping_list_path, 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        if i == 0: continue
        gene_name = line.split(';')[1]
        house_keeping_list.append(gene_name.upper() + "_1")

for record in SeqIO.parse(mouse_coding_path, "fasta"):
    gene_name: str = record.description
    sequence =  str(record.seq).upper()
    if "N" in sequence: continue
    
    class_type = PromoterTypes.CODING_HK_DNA if gene_name.upper().split()[1] in house_keeping_list else PromoterTypes.CODING_NON_HK_DNA

    sequence_set.add_sequence(
        input_dna_type=class_type,
        sequence=sequence,
    )

for record in SeqIO.parse(mouse_non_coding_path, "fasta"):
    gene_name: str = record.description
    sequence =  str(record.seq).upper()
    if "N" in sequence: continue
    sequence_set.add_sequence(
        input_dna_type=PromoterTypes.NON_CODING_NON_HK_DNA,
        sequence=sequence,
    )

with open(mouse_junk_dna_path, "r") as file:
    lines = file.readlines()
    for line in lines:
        sequence_set.add_sequence(
            input_dna_type=PromoterTypes.JUNK_DNA,
            sequence=line.strip()
        )

sequence_set.save('C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/pre_generation_sequence_set.json')
sequence_set.save_train_data('C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/processedMouseDNAData/train.h5')