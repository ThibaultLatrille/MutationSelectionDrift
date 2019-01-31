# GLOBAL IMPORTS
from codons import *
from scipy.special import binom
from scipy import integrate
import matplotlib.pyplot as plt


def normalize_bias(mut_bias):
    events = mut_bias * (2 + 1 * mut_bias)
    events += 1 * (1 + 2 * mut_bias)
    return 1.0 / events


def dfe(site_preferences, mut_bias, synonymous, subset=""):
    neighbors = generate_neighbors(synonymous, subset)

    mut_frequencies = np.power(mut_bias, nbr_weak)

    coefficients = list()
    weights = list()

    for preferences in site_preferences:
        pref_codons = np.array([preferences[codon_to_aa[i]] for i in range(len(codons))])
        codon_frequencies = mut_frequencies * pref_codons
        codon_frequencies /= np.sum(codon_frequencies)
        log_fitness = np.log(preferences)

        for i, codon_origin in enumerate(codons):
            for j, a, b in neighbors[i]:

                if synonymous:
                    coefficients.append(0.0)
                else:
                    coefficients.append(log_fitness[codon_to_aa[j]] - log_fitness[codon_to_aa[i]])

                weight = codon_frequencies[i]
                if weak_strong(nucleotides[b]) == "W":
                    weight *= mut_bias

                weights.append(weight)

    return np.array(coefficients), np.array(weights) / np.sum(weights)


def sfs(nbr_sample, site_preferences, mut_bias, theta, synonymous, subset=""):
    theta *= normalize_bias(mut_bias)
    neighbors = generate_neighbors(synonymous, subset)

    mut_frequencies = np.power(mut_bias, nbr_weak)

    sample_range = np.array(range(1, nbr_sample - 1))
    sample_sfs = np.zeros(len(sample_range))

    for preferences in site_preferences:
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

                    if synonymous:
                        pij *= 1 - x
                    else:
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

        precision = 8
        x_array = np.linspace(0, 1, 2 ** precision + 1)
        res_array = np.array([residency(x) for x in x_array])

        for index, a in enumerate(sample_range):
            y_array = res_array * np.power(x_array, a - 1) * np.power(1 - x_array, nbr_sample - a - 1)
            sample_sfs[index] += binom(nbr_sample, a) * integrate.simps(y_array, x_array)

    return sample_range, sample_sfs


def plot_dfe_sfs(alpha, at_gc_ratio, theta, sampling_effort, nbr_sites):
    bar_width = 0.35
    opacity = 1
    propensities = np.random.dirichlet(alpha * np.ones(20), nbr_sites)

    for subset in subset_list + [""]:
        param = "_".join([str(i) for i in [alpha, at_gc_ratio, theta, sampling_effort, nbr_sites, subset]])
        coeffs, weights = dfe(propensities, at_gc_ratio, synonymous=False, subset=subset)

        plt.hist(coeffs, bins=50, density=True, weights=weights,
                 label='Non-synonymous ($\\sum = {0:.3g}$)'.format(sum(weights)))

        plt.legend()
        plt.xlabel("Selection coefficient")
        plt.ylabel("Density")
        plt.title('DFE for $\\lambda = {0}$, $n = {1}$ '
                  'and $\\alpha = {2}$'.format(at_gc_ratio, sampling_effort, alpha))
        plt.savefig("../figures/dfe_{0}.png".format(param), format="png")
        plt.clf()
        plt.close('all')

        sampling, non_syn_sfs = sfs(sampling_effort, propensities, at_gc_ratio, theta, synonymous=False, subset=subset)
        plt.bar(sampling, non_syn_sfs, bar_width,
                label='Non-synonymous ($\\sum = {0:.3g}$)'.format(sum(non_syn_sfs)),
                alpha=opacity, color=BLUE)
        sampling, syn_sfs = sfs(sampling_effort, propensities, at_gc_ratio, theta, synonymous=True, subset=subset)
        plt.bar(sampling + bar_width, syn_sfs, bar_width,
                label='Synonymous ($\\sum = {0:.3g}$)'.format(sum(syn_sfs)),
                alpha=opacity, color=LIGHTGREEN)

        plt.legend()
        plt.xlabel("At derived frequency")
        plt.ylabel("Expected number of SNPs")
        plt.title('SFS for $\\lambda = {0}$, $n = {1}$ '
                  'and $\\alpha = {2}$'.format(at_gc_ratio, sampling_effort, alpha))
        plt.savefig("../figures/sfs_{0}.png".format(param), format="png")
        plt.clf()
        plt.close('all')


plot_dfe_sfs(0.1, 0.1, 0.1, 20, 100)
plot_dfe_sfs(0.1, 1.0, 0.1, 20, 100)
plot_dfe_sfs(0.1, 10., 0.1, 20, 100)
