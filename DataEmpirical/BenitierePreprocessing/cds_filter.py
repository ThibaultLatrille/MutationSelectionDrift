#!python3
import pandas as pd
from collections import Counter
import os
from Bio import SeqIO


def save_df(name, col):
    pd.DataFrame(col).to_csv(name, index=False, header=None)
    print("{0} CDS saved into '{1}'".format(len(col), name))


def coverage(filepath):
    cov = 0.0
    nbr_seqs = 0
    with open(filepath, 'r') as ali_file:
        next(ali_file)
        for line in ali_file:
            if line != "\n":
                name, seq = line.replace("  ", " ").replace("\n", "").split(" ")
                count = Counter(seq)
                cov += (count["A"] + count["C"] + count["G"] + count["T"]) / len(seq)
                nbr_seqs += 1
    return cov / nbr_seqs


for folder in ["47SP", "Vertebrates"]:
    folder_path = os.path.abspath('..') + "/" + folder
    assert (os.path.isdir(folder_path))

    for gblock in os.listdir(folder_path + '/gblocks'):
        block_path = folder_path + '/gblocks/' + gblock
        filtered_seq = [(record.name.split("_")[0], str(record.seq).upper()) for record in SeqIO.parse(block_path, "fasta")]
        ali_file = open(block_path.replace("/gblocks", "/singlegene_alignments").replace(".txt-gb", ".ali"), 'w')
        ali_file.write("{0} {1}\n".format(len(filtered_seq), len(filtered_seq[0][1])))
        ali_file.write("\n".join([" ".join(id_seq) for id_seq in filtered_seq]))
        ali_file.close()

    list_folder_cds = [i.replace(".ali", "") for i in os.listdir(folder_path + '/singlegene_alignments')]
    print("{0} CDS in the folder 'singlegene_alignments'.".format(len(list_folder_cds)))
    save_df(folder_path + "/cds.list", list_folder_cds)

    threshold = 0.99
    list_high_coverage_cds = [i for i in list_folder_cds if
                              coverage(folder_path + '/singlegene_alignments/{0}.ali'.format(i)) >= threshold]
    print("{0} CDS are with a coverage >{1}.".format(len(list_high_coverage_cds), threshold))
    save_df(folder_path + "/cds.highcoverage.list", list_high_coverage_cds)
