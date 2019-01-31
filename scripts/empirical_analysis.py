# GLOBAL IMPORTS
from analysis import *
import os
import glob

cluster = False

if cluster:
    current_dir = "/panhome/tlatrill/SimuEvol"
    cmd = "qsub"
else:
    current_dir = "/home/thibault/SimuEvol"
    cmd = "sh"

# folder = "data_empirical"
folder = "data_primates"
folder_path = "{0}/{1}".format(current_dir, folder)
fasta_files = [f for f in os.listdir(folder_path) if ".fasta" == f.strip()[-6:]]
print("Found {0} fasta files".format(len(fasta_files)))

for fasta_file in sorted(fasta_files):
    fasta_path = "{0}/{1}".format(folder_path, fasta_file)
    prefix = fasta_file.strip().replace(".fasta", "")
    print("\nAnalysis for {0}".format(prefix))

    nuc_freqs, nbr_sites, nbr_species = extract_nuc_pct(fasta_path)
    print("Number of sites: {0}".format(nbr_sites))
    print("Number of species: {0}".format(nbr_species))
    print("\t%AT of the alignment: {0:.2f}".format(nuc_freqs['A'] + nuc_freqs['T']))
    at_gc_pct_obs = (nuc_freqs['A'] + nuc_freqs['T']) / (nuc_freqs['G'] + nuc_freqs['C'])
    print("\t%AT/%GC of the alignment: {0:.2f}".format(at_gc_pct_obs))

    for hyphy_result in sorted(glob.glob("{0}/{1}*_hyout.txt".format(folder_path, prefix))):
        experiment = hyphy_result.replace(".bf_hyout.txt", "").split("_")[-1]
        hyphy_dico = dico_from_file(hyphy_result)
        hyphy_dico["n"] = nbr_sites
        format_hyphy_dico(hyphy_dico)

        print("\tResults for experiment {0}".format(experiment))
        print("\t\t w inferred: {0:.2f}".format(hyphy_dico["w"]))

        for subset in subset_list:
            omega = "w_" + subset
            if omega in hyphy_dico:
                print("\t\t {0} inferred: {1:.2f}".format(omega, hyphy_dico[omega]))

        if ("pnG" in hyphy_dico) and ("pnG" in hyphy_dico):
            gc_pct = hyphy_dico["pnG"] + hyphy_dico["pnC"]
            print("\t\t %AT inferred from the mutational process: {0:.2f}".format(1 - gc_pct))
            print("\t\t %AT/%GC inferred from the mutational process: {0:.2f}".format((1 - gc_pct) / gc_pct))

        print("\t\t %AT inferred from the Mut-Sel equilibrium: {0:.2f}".format(equilibrium_at_pct(hyphy_dico)))
        print("\t\t %AT/%GC inferred from the Mut-Sel equilibrium: {0:.2f}".format(equilibrium_lambda(hyphy_dico)))
