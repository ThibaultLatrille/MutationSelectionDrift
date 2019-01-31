# GLOBAL IMPORTS
from codons import *


def dico_from_file(filename):
    tmp_dico = {}
    tmp_file = open(filename, "r")
    for line in tmp_file:
        split_line = line.split("=")
        if len(split_line) > 1:
            value = split_line[1].strip()
            try:
                tmp_dico[split_line[0]] = float(value)
            except:
                pass
    tmp_file.close()
    return tmp_dico


def extract_nuc_pct(fasta_path):
    ''' Compute nucleotides frequencies from fasta file '''
    nucindex = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
    total = 0.0
    nbr_sites = 0
    nbr_species = 0
    fasta_file = open(fasta_path, 'r')
    for seq_unstriped in fasta_file:
        if seq_unstriped[0] != ">":
            nbr_species += 1
            seq = seq_unstriped.strip()
            assert len(seq) % 3 == 0
            nbr_sites = int(len(seq) / 3)
            for site in seq:
                if site in nucleotides:
                    nucindex[site] += 1.
                    total += 1.
    fasta_file.close()

    return {k: v / total for k, v in nucindex.items()}, nbr_sites, nbr_species


def format_hyphy_dico(hyphy_dico):
    if "pnCG" in hyphy_dico:
        hyphy_dico["pnC"] = hyphy_dico["pnCG"]
        hyphy_dico["pnG"] = hyphy_dico["pnCG"]
        hyphy_dico["pnA"] = 0.5 - hyphy_dico["pnCG"]
    hyphy_dico["pnT"] = 1.0 - (hyphy_dico["pnA"] + hyphy_dico["pnC"] + hyphy_dico["pnG"])

    if "epsA" in hyphy_dico and "epsM" not in hyphy_dico:
        hyphy_dico["epsM"] = 20 - sum([hyphy_dico["eps" + aa] for aa in amino_acids_set if aa != "M"])

    if 'w' not in hyphy_dico:
        codon_frequencies = np.ones(len(codons))
        for codon_index, codon in enumerate(codons):
            for nuc in codon:
                codon_frequencies[codon_index] *= hyphy_dico["pn" + nuc]
            epsilon = "eps" + codon_table[codon]
            if epsilon in hyphy_dico:
                codon_frequencies[codon_index] *= hyphy_dico[epsilon]

        codon_frequencies /= np.sum(codon_frequencies)

        d_dict, d0_dict = defaultdict(float), defaultdict(float)
        d, d0 = 0.0, 0.0

        for x, codon_origin in enumerate(codons):
            for y, a, b in non_syn_neighbors[x]:
                nuc_origin, nuc_target = nucleotides[a], nucleotides[b]
                codon_target = codons[y]

                mut_flow_tmp = codon_frequencies[x] * hyphy_dico["pn" + nuc_target]

                exchan = "exch" + "".join(sorted((nuc_target, nuc_origin)))
                if exchan in hyphy_dico:
                    mut_flow_tmp *= hyphy_dico[exchan]

                p_fix = 1.0
                beta = "b_" + "".join(sorted((codon_table[codon_origin], codon_table[codon_target])))
                if beta in hyphy_dico:
                    if hyphy_dico[beta] > 100.0:
                        print("{0}={1} ({2} sites)".format(beta, hyphy_dico[beta], hyphy_dico["n"]))
                        hyphy_dico[beta] = 1.0
                    p_fix *= hyphy_dico[beta]

                omega_subset = "w_" + weak_strong(nuc_origin) + weak_strong(nuc_target)
                if omega_subset in hyphy_dico:
                    p_fix *= hyphy_dico[omega_subset]

                epsilon = "eps" + codon_table[codon_target]
                if epsilon in hyphy_dico:
                    p_fix *= hyphy_dico[epsilon]

                d += mut_flow_tmp * p_fix
                d0 += mut_flow_tmp

                if omega_subset not in hyphy_dico:
                    d_dict[omega_subset] += mut_flow_tmp * p_fix
                    d0_dict[omega_subset] += mut_flow_tmp

        hyphy_dico["w"] = d / d0
        for omega_subset in d_dict.keys():
            if d_dict[omega_subset] != 0.0 and (omega_subset not in hyphy_dico):
                hyphy_dico[omega_subset] = d_dict[omega_subset] / d0_dict[omega_subset]


def equilibrium_lambda(hyphy_dico):
    at_pct = equilibrium_at_pct(hyphy_dico)
    return at_pct / (1 - at_pct)


def equilibrium_at_pct(hyphy_dico):
    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    codon_frequencies = np.ones(len(codons))

    for index, codon in enumerate(codons):
        for nuc in codon:
            codon_frequencies[index] *= hyphy_dico["pn" + nuc]
        epsilon = "eps" + codon_table[codon]
        if epsilon in hyphy_dico:
            codon_frequencies[index] *= hyphy_dico[epsilon]

    codon_frequencies /= np.sum(codon_frequencies)

    return np.sum(codon_frequencies * nbr_weak) / 3


def prefs_path_to_list(preferences_path):
    preferences_list = []
    preferences_file = open(preferences_path, "r")
    preferences_file.readline()
    for line in preferences_file:
        preferences = list(map(float, line.strip().split(" ")[3:]))
        assert len(preferences) == 20
        preferences_list.append(preferences)
    preferences_file.close()
    return preferences_list


def theoretical_at_gc_pct(preferences_list, mut_bias):
    at_pct = theoretical_at_pct(preferences_list, mut_bias)
    return at_pct / (1 - at_pct)


def theoretical_at_pct(preferences_list, mut_bias):
    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    at_pct_list = []

    for preferences in preferences_list:
        codon_frequencies = np.power(mut_bias, nbr_weak)
        pref_codons = [preferences[codon_to_aa[i]] for i in range(len(codons))]
        codon_frequencies *= np.power(pref_codons, 1)
        codon_frequencies /= np.sum(codon_frequencies)

        at_pct_list.append(np.sum(codon_frequencies * nbr_weak) / 3)

    return np.mean(at_pct_list)