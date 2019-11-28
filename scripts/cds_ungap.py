#!python3
import os
from Bio import AlignIO, SeqIO, SeqRecord, Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Data.CodonTable import TranslationError
from ete3 import Tree

f = os.path.abspath('..') + "/DataEmpirical/Cetacea"

phy = AlignIO.read("{0}/datadryad/DATASET_B.phylip".format(f), format="phylip-relaxed", alphabet=IUPACAmbiguousDNA())
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

        seqs = []
        for rec in phy[:, down - 1:up]:
            if rec.id not in taxa:
                continue

            seq = Seq.Seq(str(rec.seq).replace("-", ""))
            if len(seq) % 3 != 0 or len(seq) == 0:
                continue

            try:
                seq[:-3].translate()
            except TranslationError:
                continue

            seqs.append(SeqRecord.SeqRecord(seq, rec.id, "", ""))

        if len(seqs) < len(phy) / 3:
            continue
        SeqIO.write(seqs, "{0}/fasta/{1}.fasta".format(f, name), 'fasta')