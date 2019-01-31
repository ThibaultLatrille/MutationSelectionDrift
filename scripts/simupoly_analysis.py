from glob import glob
from csv import reader
from math import isfinite
from numpy import nanmean, nanstd
from codons import *
from scipy.special import binom
from scipy import integrate
import matplotlib.pyplot as plt

folder = "data_sfs"
folder_path = "/home/thibault/SimuEvol/{0}".format(folder)


def sfs_synonymous(nbr_sample, nbr_sites, site_mutation_rate, mut_bias, population_size):
    theta = 4 * population_size * nbr_sites * site_mutation_rate
    neighbors = generate_neighbors(True)

    mut_frequencies = np.power(mut_bias, nbr_weak)
    mut_frequencies /= np.sum(mut_frequencies)

    sample_sfs = np.zeros(nbr_sample + 1)

    def residency(x):
        res = 0
        for i, codon_origin in enumerate(codons):
            tmp_res = 0
            for j, a, b in neighbors[i]:
                pij = theta

                pij *= 1 - x

                if weak_strong(nucleotides[b]) == "W":
                    pij *= mut_bias

                tmp_res += pij

            res += mut_frequencies[i] * tmp_res
        return res

    precision = 6
    x_array = np.linspace(0, 1, 2 ** precision + 1)
    res_array = np.array([residency(x) for x in x_array])

    for a in range(1, nbr_sample):
        y_array = res_array * np.power(x_array, a - 1) * np.power(1 - x_array, nbr_sample - a - 1)
        sample_sfs[a] += binom(nbr_sample, a) * integrate.simps(y_array, x_array)

    return sample_sfs


def sfs_non_synonymous(nbr_sample, site_preferences, site_mutation_rate, mut_bias, population_size):
    theta = 4 * population_size * site_mutation_rate
    neighbors = generate_neighbors(False)

    mut_frequencies = np.power(mut_bias, nbr_weak)

    sample_sfs = np.zeros(nbr_sample + 1)

    for n, preferences in enumerate(site_preferences):
        print("{0:.2f}% computed".format(100 * n / len(site_preferences)))
        pref_codons = np.array([preferences[codon_to_aa[i]] for i in range(len(codons))])
        codon_frequencies = mut_frequencies * pref_codons
        codon_frequencies /= np.sum(codon_frequencies)
        log_fitness = np.log(preferences)

        def residency(x):
            res = 0
            for i, codon_origin in enumerate(codons):
                tmp_res = 0
                for j, a, b in neighbors[i]:
                    pij = theta

                    s = log_fitness[codon_to_aa[j]] - log_fitness[codon_to_aa[i]]
                    if s == 0:
                        pij *= 1 - x
                    else:
                        pij *= (1 - np.exp(- s * (1 - x))) / (1 - np.exp(-s))

                    if weak_strong(nucleotides[b]) == "W":
                        pij *= mut_bias

                    tmp_res += pij

                res += codon_frequencies[i] * tmp_res
            return res

        precision = 6
        x_array = np.linspace(0, 1, 2 ** precision + 1)
        res_array = np.array([residency(x) for x in x_array])

        for a in range(1, nbr_sample):
            y_array = res_array * np.power(x_array, a - 1) * np.power(1 - x_array, nbr_sample - a - 1)
            sample_sfs[a] += binom(nbr_sample, a) * integrate.simps(y_array, x_array)

    return sample_sfs


for tsv_path in sorted(glob("{0}/np*.tsv".format(folder_path))):
    with open(tsv_path, 'r') as tsvfile:
        tsvin = list(reader(tsvfile, delimiter='\t'))

        n = tsv_path.split("/")[-1].replace('.tsv', '').split("_")
        protein = n[0]
        mu = float(n[1])
        population_size = int(n[2])
        mut_bias = float(n[3])
        beta = float(n[4])
        title = "{0} $\\mu={1}$, $Ne={2}$, $\\lambda={3}$, $\\beta={4}$".format(protein, mu, population_size, mut_bias,
                                                                                beta)

        if protein in ['gal4', 'lactamase', 'np', 'ha']:
            import os

            preferences_list = []
            file_path = '../data_prefs/{0}.txt'.format(protein)
            if os.path.isfile(file_path):
                with open(file_path, 'r') as prefs_file:
                    prefs_file.readline()
                    for line in prefs_file:
                        prefs = [float(f) for f in line.replace('\n', '').replace('\t', ' ').split(' ')[3:]]
                        assert len(prefs) == 20
                        preferences_list.append(np.array(prefs) / np.sum(prefs))
            propensities = np.power(np.array(preferences_list), beta)
        else:
            propensities = np.random.dirichlet(np.ones(len(amino_acids)), int(protein))

        print("{0} sites".format(propensities.shape[0]))

        histogram = [row for row in tsvin if len(row) > 3 and row[2] == "0"]
        dict_label_histogram = {}
        for row in histogram:
            if not row[3] in dict_label_histogram:
                dict_label_histogram[row[3]] = []
            histo = [float(i) for i in row[4].split(" ")]
            dict_label_histogram[row[3]].append(histo)
        bar_width = 0.35
        opacity = 1
        for label, hist in dict_label_histogram.items():
            assert len(set([len(r) for r in hist])) == 1
            my_dpi = 92
            plt.figure(figsize=(800 / my_dpi, 600 / my_dpi), dpi=my_dpi)
            sfs_hist = np.mean(np.array(hist), axis=0)
            nbr_sample = len(sfs_hist) - 1
            range_sample = np.array(range(nbr_sample + 1))
            print("Computing theoretical SFS")
            if label == "SFSs":
                predicted = sfs_synonymous(nbr_sample, len(propensities), mu, mut_bias, population_size)
            else:
                predicted = sfs_non_synonymous(nbr_sample, propensities, mu, mut_bias, population_size)
            plt.bar(range_sample, predicted, bar_width,
                    label='Predicted {0} ($\\sum = {1:.3g}$)'.format(label, sum(predicted)),
                    alpha=opacity, color=RED)
            plt.bar(range_sample + bar_width, sfs_hist, bar_width,
                    label='{0} ($\\sum = {1:.3g}$)'.format(label, sum(sfs_hist)),
                    alpha=opacity, color=BLUE)
            plt.legend()
            plt.yscale("log")
            plt.title(title)
            plt.tight_layout()
            plt.savefig("../figures/{0}_{1}".format(tsv_path.split("/")[-1].replace("tsv", "png"), label), format="png")
            plt.show()
            plt.clf()
            plt.close('all')

        data = [row for row in tsvin if len(row) > 3 and row[2] != "0"]
        dict_num_label_data = {}
        for row in data:
            if not row[2] in dict_num_label_data:
                dict_num_label_data[row[2]] = {}
            if not row[3] in dict_num_label_data[row[2]]:
                dict_num_label_data[row[2]][row[3]] = {"x": [], "y": []}

            if isfinite(float(row[4])):
                dict_num_label_data[row[2]][row[3]]["x"].append(float(row[1]))
                dict_num_label_data[row[2]][row[3]]["y"].append(float(row[4]))

        for num, dict_label_data in dict_num_label_data.items():
            my_dpi = 92
            plt.figure(figsize=(800 / my_dpi, 600 / my_dpi), dpi=my_dpi)
            for label, dict_data in dict_label_data.items():
                plt.plot(dict_data["x"], dict_data["y"], label=label)
                print("{0}={1:.3f}±{2:.3f}".format(label, nanmean(dict_data["y"]), nanstd(dict_data["y"])))
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.xlabel('$t$', fontsize=24)
            plt.legend()
            plt.title(title)
            plt.tight_layout()
            plt.savefig("../figures/{0}_{1}".format(tsv_path.split("/")[-1].replace("tsv", "png"), num),
                        format="png")
            plt.show()
            plt.clf()
            plt.close('all')
