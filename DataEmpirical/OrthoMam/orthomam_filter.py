#!python3
import pandas as pd
from collections import Counter
import os


def save_df(name, col):
    col.to_csv("./" + name, index=False, header=None)
    print("{0} CDS saved into '{1}'".format(len(col), name))


df_orthomam_cds = pd.read_csv("./orthomam.tsv", sep="\t")
col_orthomam_cds = df_orthomam_cds["gene"]
print("{0} CDS with mutsel inference.".format(len(col_orthomam_cds)))

col_pursel_cds = df_orthomam_cds.query("pp_m2a < 0.1 and pp_mutsel < 0.1 and pp_mutselfix < 0.1")["gene"]
print("{0} CDS with mutsel inference not under positive selection".format(len(col_pursel_cds)))

list_folder_cds = [i.replace(".ali", "") for i in os.listdir('./singlegene_alignments')]
print("{0} CDS in the folder 'singlegene_alignments'.".format(len(list_folder_cds)))

col_filtered_cds = col_pursel_cds[col_pursel_cds.isin(list_folder_cds)]
print("{0} CDS in the folder 'singlegene_alignments' and not under positive selection.".format(len(col_filtered_cds)))
save_df("cds.filtered.list", col_filtered_cds)


def coverage(filepath):
    cov = 0.0
    nbr_seqs = 0
    with open(filepath, 'r') as ali_file:
        next(ali_file)
        for line in ali_file:
            if line != "\n":
                name, seq = line.replace("  ", " ").replace("\n", "").split(" ")
                count = Counter(seq.upper())
                cov += (count["A"] + count["C"] + count["G"] + count["T"]) / len(seq)
                nbr_seqs += 1
    return cov / nbr_seqs


threshold = 0.99
list_high_coverage_cds = [i for i in list_folder_cds if coverage('./singlegene_alignments/{0}.ali'.format(i)) >= threshold]
col_filtered_high_coverage_cds = col_filtered_cds[col_filtered_cds.isin(list_high_coverage_cds)]
print("{0} CDS are not under positive selection and with a coverage >{1}.".format(len(col_filtered_high_coverage_cds), threshold))
save_df("cds.filtered.highcoverage.list", col_filtered_high_coverage_cds)

col_high_coverage_cds = col_orthomam_cds[col_orthomam_cds.isin(list_high_coverage_cds)]
print("{0} CDS are with a coverage >{1}.".format(len(col_high_coverage_cds), threshold))
save_df("cds.highcoverage.list", col_high_coverage_cds)
