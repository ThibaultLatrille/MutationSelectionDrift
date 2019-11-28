#!python3
import pandas as pd
from collections import Counter
import os
from Bio import SeqIO


def save_df(name, col):
    pd.DataFrame(col).to_csv(name, index=False, header=None)
    print("{0} CDS saved into '{1}'".format(len(col), name))


def size(filepath):
    s = 0.0
    with open(filepath, 'r') as ali_read:
        s = int(ali_read.readline().split(" ")[1])
        assert (s > 0)
    return s


def coverage(filepath):
    cov = 0.0
    nbr_seqs = 0
    with open(filepath, 'r') as ali_read:
        next(ali_read)
        for line in ali_read:
            if line != "\n":
                name, seq = line.replace("  ", " ").replace("\n", "").split(" ")
                count = Counter(seq)
                cov += (count["A"] + count["C"] + count["G"] + count["T"]) / len(seq)
                nbr_seqs += 1
    return cov / nbr_seqs


aln_folder = "gblocks"
folder_path = os.path.abspath('..') + "/DataEmpirical/Drosophila_v2"
assert (os.path.isdir(folder_path))

for aln in os.listdir(folder_path + '/' + aln_folder):
    block_path = folder_path + '/' + aln_folder + '/' + aln
    filtered_seq = [(record.name.split("_")[0], str(record.seq).upper().replace("!", "-")) for record in
                    SeqIO.parse(block_path, "fasta")]
    seq_sizes = set([len(v) for k, v in filtered_seq])
    assert (len(seq_sizes) == 1)
    if seq_sizes.pop() % 3 != 0:
        print(aln + " is not a multiple of 3")
        continue
    ali_file = open(block_path.replace("/" + aln_folder, "/singlegene_alignments") + ".ali", 'w')
    ali_file.write("{0} {1}\n".format(len(filtered_seq), len(filtered_seq[0][1])))
    ali_file.write("\n".join([" ".join(id_seq) for id_seq in filtered_seq]))
    ali_file.close()

list_folder_cds = [i.replace(".ali", "") for i in os.listdir(folder_path + '/singlegene_alignments')]
print("{0} CDS in the folder 'singlegene_alignments'.".format(len(list_folder_cds)))
save_df(folder_path + "/cds.list", list_folder_cds)

for size_min in [300, 600, 900, 1200]:
    size_max = size_min + 300
    list_long_cds = [i for i in list_folder_cds if
                     size_max >= size(folder_path + '/singlegene_alignments/{0}.ali'.format(i)) >= size_min]
    print("{0} CDS are at least {1}bp long.".format(len(list_long_cds), size_min))
    save_df(folder_path + "/cds.size{0}.list".format(size_min), list_long_cds)

for threshold in []:
    list_high_coverage_cds = [i for i in list_folder_cds if
                              coverage(folder_path + '/singlegene_alignments/{0}.ali'.format(i)) >= threshold]
    print("{0} CDS are with a coverage >{1}.".format(len(list_high_coverage_cds), threshold))
    save_df(folder_path + "/cds.coverage{0}.list".format("{0:.2f}".format(threshold).split('.')[1]),
            list_high_coverage_cds)
