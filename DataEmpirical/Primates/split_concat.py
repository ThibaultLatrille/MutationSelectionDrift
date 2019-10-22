import os

ali_path = "./CDS.ali"
os.makedirs("./singlegene_alignments", exist_ok=True)

ali_dict = dict()
with open(ali_path, 'r') as ali_file:
    next(ali_file)
    for line in ali_file:
        if line != "\n":
            name, seq = line.replace("\t", " ").replace("  ", " ").replace("\n", "").split(" ")
            ali_dict[name] = seq


seq_sizes = set([len(v) for v in ali_dict.values()])
assert (len(seq_sizes) == 1)
seq_size = seq_sizes.pop()
exon = 1
exon_size = 900
assert(exon_size % 3 == 0)
start_site = 0
while start_site < seq_size:
    end_site = min(start_site + exon_size, seq_size)
    with open("./singlegene_alignments/CDS_{0}.ali".format(exon), 'w') as ali_file:
        ali_file.write("{0} {1}\n".format(len(ali_dict), end_site - start_site))
        ali_file.write("\n".join(["{0} {1}".format(k, v[start_site:end_site]) for k, v in ali_dict.items()]))
    start_site += exon_size
    exon += 1

