#!python3
import pandas as pd
import numpy as np
import os
from Bio import AlignIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Data.CodonTable import TranslationError
from ete3 import Tree

f = os.path.abspath('..') + "/DataEmpirical/Cetacea"

phy = AlignIO.read("{0}/datadryad/DATASET_B.phylip".format(f), format="phylip-relaxed", alphabet=IUPACUnambiguousDNA())
print("{0} taxa.".format(len(phy)))
taxa = Tree("{0}/rootedtree.nhx".format(f), format=1).get_leaf_names()

precision_dict = {}
coverage_dict = {}
with open("{0}/datadryad/Cetacea_gene_partition.txt".format(f), "r") as gene_partition:
    for line in gene_partition:
        name, pos = line.replace("DNA,", "").replace(" ", "").split("=")
        down, up = pos.split("-")
        down, up = int(down), int(up)
        diff = 1 + up - down
        if diff % 3 != 0:
            continue

        sequences = phy[:, down - 1:up]
        output = phy[:, :0]
        filtered = [rec for rec in sequences if rec.id in taxa]
        for pos in range(0, int(diff / 3)):
            keep_site = True
            for sr in filtered:
                site = sr.seq[pos * 3:(pos + 1) * 3]
                try:
                    site.translate(gap="-")
                except TranslationError:
                    try:
                        site.translate(gap="?")
                    except TranslationError:
                        keep_site = False
                        break
            if keep_site:
                output += phy[:, pos * 3:(pos + 1) * 3]
        seq_size = len(output[0])
        if seq_size < 10:
            continue

        precision_dict[name] = seq_size / diff
        ids_seqs = {fasta.id: str(fasta.seq) for fasta in output if (fasta.id in taxa)}
        seqs_cov = {k: (len(v) - v.count("-") - v.count("?")) / len(v) for k, v in ids_seqs.items()}
        for k, v in seqs_cov.items():
            if v == 0:
                print("Removing " + k)
                ids_seqs.pop(k)
        coverage_dict[name] = sum(seqs_cov.values()) / len(ids_seqs)
        assert(len(set([len(s) for s in ids_seqs.values()])) == 1)
        ali_file = open("{0}/singlegene_alignments/{1}.ali".format(f, name), 'w')
        ali_file.write("{0} {1}\n".format(len(ids_seqs), seq_size))
        ali_file.write("\n".join([" ".join(k_v) for k_v in ids_seqs.items()]))
        ali_file.close()
        print("{0}: {1:.1f}% acc; {2:.1f}% cov".format(name, precision_dict[name] * 100, coverage_dict[name] * 100))

for precision in np.arange(0.9, 1.0, 0.02):
    for coverage in np.arange(0.7, 1.0, 0.05):
        cds = [k for k, v in precision_dict.items() if v >= precision and coverage_dict[k] > coverage]
        print("{0} CDS are with a acc > {1} and coverage > {2}.".format(len(cds), precision, coverage))
        filename = f + "/cds.{0}acc.{1}cov.list".format("{0:.2f}".format(precision).split('.')[1],
                                                        "{0:.2f}".format(coverage).split('.')[1])
        pd.DataFrame(cds).to_csv(filename, index=False, header=None)
        print("{0} CDS saved into '{1}'".format(len(cds), filename))
