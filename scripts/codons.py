# GLOBAL IMPORTS
import numpy as np
from collections import defaultdict

subset_list = ["WS", "WW", "SS", "SW"]
codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}

nucleotides = "ACGT"
weak_nucleotides = "AT"
strong_nucleotides = "CG"
nucindex = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

assert len(codon_table.keys()) == 64, "There is not 64 codons in the codon table"
codons = sorted(list([k for k, v in codon_table.items() if v != 'X']))
assert len(codons) == 61, "There is not 3 stop codons in the codon table"
assert not [_nuc for _nuc in "".join(codons) if
            (_nuc not in nucleotides)], "There is a codon with an unrecognized nucleotide"
nbr_weak = np.array([_c.count("A") + _c.count("T") for _c in codons])
nbr_strong = 3 - nbr_weak
position_nbr_weak = dict()
for _p in range(1, 4):
    position_nbr_weak[_p] = np.array([_c[_p - 1].count("A") + _c[_p - 1].count("T") for _c in codons])

amino_acids_set = set(codon_table.values())
amino_acids_set.remove('X')

aa_char_to_int = {v: k for k, v in enumerate(sorted(amino_acids_set))}
amino_acids = "".join(amino_acids_set)
assert len(amino_acids) == 20, "There is not 20 amino-acids in the codon table"


def generate_aa_table():
    aa_dict = {}
    for index, codon in enumerate(codons):
        aa = codon_table[codon]
        if aa not in aa_dict:
            aa_dict[aa] = []
        aa_dict[aa].append(codon)
    assert len(aa_dict) == 20, "There is not 20 amino-acids in the aa table"
    return aa_dict


def weak_strong(nuc):
    if nuc == "A" or nuc == "T":
        return "W"
    else:
        return "S"


def generate_c_to_aa():
    codon_to_aa_list = [None] * len(codons)
    for codon_index, codon in enumerate(codons):
        codon_to_aa_list[codon_index] = aa_char_to_int[codon_table[codon]]
    return codon_to_aa_list


def generate_neighbors(synonymous=True, subset=""):
    neighbor_dict = dict()
    for x, codon_origin in enumerate(codons):
        neighbor_dict[x] = []
        for y, codon_target in enumerate(codons):
            mutations = [(n_from, n_to) for n_from, n_to in zip(codon_origin, codon_target) if n_from != n_to]
            if len(mutations) == 1:
                if synonymous == (codon_table[codon_origin] == codon_table[codon_target]):
                    nuc_origin, nuc_target = mutations[0]
                    if (subset == "") or ((weak_strong(nuc_origin) + weak_strong(nuc_target)) == subset):
                        neighbor_dict[x].append((y, nucleotides.index(nuc_origin), nucleotides.index(nuc_target)))
    return neighbor_dict


aa_table = generate_aa_table()
codon_to_aa = generate_c_to_aa()
non_syn_neighbors = generate_neighbors(synonymous=False)
syn_neighbors = generate_neighbors(synonymous=True)
all_neighbors = {k: non_syn_neighbors[k] + syn_neighbors[k] for k in syn_neighbors.keys()}
assert len(all_neighbors) == 61


def nested_dict_init():
    return defaultdict(nested_dict_init)


RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"


def get_color(i):
    colors = [RED, YELLOW, BLUE, GREEN, LIGHTGREEN]
    return colors[i % len(colors)]
