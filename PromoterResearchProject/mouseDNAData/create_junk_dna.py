from Bio import SeqIO

chX_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/mouseDNAData/chrX.fa'
save_path = 'C:/Users/kunal/Documents/BoekeLabResearch/PromoterResearchProject/mouseDNAData/junk_dna.txt'

for record in SeqIO.parse(chX_path, "fasta"):
    print(f"ID: {record.id}")
    print(f"Description: {record.description}")
    junk_region =  record.seq[11325992:11954792]
    print("N" in junk_region)

junk_region = str(junk_region).upper()

sequence_size = 1024
save_string = ""
for i in range(len(junk_region)//1024):
    save_string += junk_region[sequence_size * i: sequence_size * (i + 1)] + "\n"

with open(save_path, "w") as file:
    file.write(save_string)