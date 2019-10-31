import os

base_path = "/Users/FlaviuVadan/Documents/flav/current_courses/BINF400/proteoglycan-pathway-evolution/src/data"
gene_path = os.path.join(base_path, "organisms/Acan_ENSG00000157766.txt")
organisms = []
with open(gene_path, 'r') as f:
    for l in f.readlines():
        if l.startswith(">"):
            organisms.append(l.replace(">", "").strip())

for org in organisms:
    os.system("bash parser.sh {}".format(org))
    break

with open(os.path.join(base_path, "alignments/Acan_msa.temp"), "r") as f:
    alignment = f.read()
    print(alignment.replace("homo_sapiens", "").strip().replace(" ", "").replace("\n", ""))
