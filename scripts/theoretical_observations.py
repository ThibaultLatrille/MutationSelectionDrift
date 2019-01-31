# GLOBAL IMPORTS
from codons import *
import matplotlib.pyplot as plt


def simpson_diversity(frequencies):
    return 1.0 / np.sum(frequencies * frequencies)


def shanon_diversity(frequencies):
    return np.exp(-np.sum(frequencies * np.log(frequencies)))


def hill_diversity(frequencies, q):
    if q == 1:
        return shanon_diversity(frequencies)
    else:
        return np.power(np.sum(np.power(frequencies, q)), 1.0 / (1.0 - q))


def compute_diversity(frequencies_seq):
    n = frequencies_seq.shape[0]
    diversity_seq = np.zeros(n)
    for site in range(n):
        diversity_seq[site] = simpson_diversity(frequencies_seq[site, :])

    diversity_mean = simpson_diversity(np.mean(frequencies_seq, axis=0))
    return np.mean(diversity_seq), diversity_mean


def compute_diversity_aa(codon_frequencies_seq):
    n = codon_frequencies_seq.shape[0]
    aa_frequencies_seq = np.zeros((n, len(amino_acids)))
    for site in range(n):
        for codon_index, freq in enumerate(codon_frequencies_seq[site, :]):
            aa_frequencies_seq[site, codon_to_aa[codon_index]] += freq

    return compute_diversity(aa_frequencies_seq)


def compute_lambda_obs(codon_frequencies_seq, position=0):
    n = codon_frequencies_seq.shape[0]
    at_pct_seq = np.zeros(n)
    for site in range(n):
        codon_frequencies_site = codon_frequencies_seq[site, :]
        if position == 0:
            at_pct_seq[site] = np.sum(codon_frequencies_site * nbr_weak) / 3
        else:
            at_pct_seq[site] = np.sum(codon_frequencies_site * position_nbr_weak[position])

    at_pct = np.mean(at_pct_seq)
    return at_pct / (1 - at_pct)


def compute_omega(codon_frequencies_seq, aa_fitness_seq, mutation_matrix, gBGC_param):
    sub_flow_non_syn_subset, mut_flow_non_syn_subset = defaultdict(float), defaultdict(float)
    omega_subset = defaultdict(float)

    n = codon_frequencies_seq.shape[0]
    for site in range(n):
        aa_fitness_site = aa_fitness_seq[site, :]
        codon_frequencies_site = codon_frequencies_seq[site, :]
        for x in range(len(codons)):

            for y, a, b in non_syn_neighbors[x]:
                sel_coef = aa_fitness_site[codon_to_aa[y]] - aa_fitness_site[codon_to_aa[x]]

                if nucleotides[b] in strong_nucleotides:
                    sel_coef += gBGC_param

                if abs(sel_coef) < 1e-10:
                    p_fix = 1.0
                else:
                    p_fix = sel_coef / (1. - np.exp(-sel_coef))

                mut_flow_tmp = codon_frequencies_site[x] * mutation_matrix[a][b]

                subset = weak_strong(nucleotides[a]) + weak_strong(nucleotides[b])
                sub_flow_non_syn_subset[subset] += mut_flow_tmp * p_fix
                mut_flow_non_syn_subset[subset] += mut_flow_tmp

                sub_flow_non_syn_subset["all"] += mut_flow_tmp * p_fix
                mut_flow_non_syn_subset["all"] += mut_flow_tmp

    for subset in mut_flow_non_syn_subset.keys():
        omega_subset[subset] = (sub_flow_non_syn_subset[subset] / mut_flow_non_syn_subset[subset])

    return omega_subset


def generate_frequencies(beta_param, lambda_param, gBGC_param, propensities_seq):
    mutation_matrix = np.ones((len(nucleotides), len(nucleotides)))
    for a in range(len(nucleotides)):
        for b in range(len(nucleotides)):
            if nucleotides[b] in weak_nucleotides:
                mutation_matrix[a, b] = lambda_param

    aa_propensities_seq = np.power(propensities_seq, beta_param)
    n, _ = aa_propensities_seq.shape
    aa_fitness_seq = np.log(aa_propensities_seq)
    codon_frequencies_seq = np.zeros((n, len(codons)))
    mut_frequencies = np.power(lambda_param, nbr_weak)
    mut_frequencies *= np.exp(gBGC_param * nbr_strong)

    for site in range(n):
        aa_propensities_site = aa_propensities_seq[site, :]
        codon_propensities_site = np.array([aa_propensities_site[codon_to_aa[i]] for i in range(len(codons))])
        codon_frequencies_site = mut_frequencies * codon_propensities_site
        codon_frequencies_seq[site, :] = codon_frequencies_site / np.sum(codon_frequencies_site)

    omega = compute_omega(codon_frequencies_seq, aa_fitness_seq, mutation_matrix, gBGC_param)
    lambda_dict = dict()
    for position in range(4):
        lambda_dict[position] = compute_lambda_obs(codon_frequencies_seq, position)
    diversity_codon, diversity_mean_codon = compute_diversity(codon_frequencies_seq)
    diversity_aa, diversity_mean_aa = compute_diversity_aa(codon_frequencies_seq)

    return omega["all"], omega["SS"], omega["SW"], omega["WS"], omega["WW"], \
           lambda_dict[0], lambda_dict[1], lambda_dict[2], lambda_dict[3], \
           diversity_codon, diversity_mean_codon, diversity_aa, diversity_mean_aa


def plot_summary_stats(lambda_axis, beta_axis, b_axis, protein, dimensions=list([1])):
    params_axis_list = [lambda_axis, beta_axis, b_axis]
    params_axis_len_list = [len(x) for x in params_axis_list]
    labels = ["\\lambda", "N_{e}", "B"]
    assert min(params_axis_len_list) == 1, "One of the parameter axis must have length 1"
    assert len(set(params_axis_len_list)) == 3, "The parameter axis must have different lengths"
    order_to_display = {}
    x_axis, y_axis = [], []
    z_param = 0
    for param_index, param_axis in enumerate(params_axis_list):
        if len(param_axis) == max(params_axis_len_list):
            x_axis = param_axis
            order_to_display[param_index] = 0
        elif len(param_axis) != 1:
            y_axis = param_axis
            order_to_display[param_index] = 1
        else:
            z_param = param_axis[0]
            order_to_display[param_index] = 2

    display_to_order = {v: k for k, v in order_to_display.items()}
    assert len(x_axis) > 1 and len(y_axis) > 1

    stat_axis = ['omega', 'omega_SS', 'omega_SW', 'omega_WS', 'omega_WW',
                 'at_over_gc', 'at_over_gc_1', 'at_over_gc_2', 'at_over_gc_3',
                 'diversity_site_codon', 'diversity_seq_codon',
                 'diversity_site_aa', 'diversity_seq_aa']

    stat_names = ['$\\omega$',
                  '$\\omega_{\\mathrm{GC} \\rightarrow \\mathrm{GC}}$',
                  '$\\omega_{\\mathrm{GC} \\rightarrow \\mathrm{AT}}$',
                  '$\\omega_{\\mathrm{AT} \\rightarrow \\mathrm{GC}}$',
                  '$\\omega_{\\mathrm{AT} \\rightarrow \\mathrm{AT}}$',
                  'AT/GC', 'AT/GC', 'AT/GC', 'AT/GC',
                  'Diversity', 'Diversity',
                  'Diversity', 'Diversity']

    array_dict = dict()
    for stat_index in range(len(stat_axis)):
        array_dict[stat_index] = np.zeros((len(y_axis), len(x_axis)))

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
        propensities = np.array(preferences_list)
    else:
        propensities = np.random.dirichlet(np.ones(len(amino_acids)), int(protein))
    print("{0} sites".format(propensities.shape[0]))

    for y_index, y_param in enumerate(y_axis):
        for x_index, x_param in enumerate(x_axis):
            display_list = [x_param, y_param, z_param]
            lambda_param, beta_param, b_param = display_list[order_to_display[0]], \
                                                display_list[order_to_display[1]], \
                                                display_list[order_to_display[2]]
            stat_list = generate_frequencies(beta_param, lambda_param, b_param, propensities)
            for stat_index, stat in enumerate(stat_list):
                array_dict[stat_index][y_index][x_index] = stat

            pct = (y_index * len(x_axis) + (x_index + 1)) / (len(y_axis) * len(x_axis))
            print("{0:.4g}% computed".format(100 * pct))

    for stat_index, x_y_array in array_dict.items():

        for dimension in dimensions:
            my_dpi = 196
            plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
            if 1 == dimension:
                for y_index, z_axis in enumerate(x_y_array):
                    plt.plot(x_axis, z_axis, color=get_color(y_index),
                             linewidth=2,
                             label="${0}={1:.2f}$".format(labels[display_to_order[1]], y_axis[y_index]))
                if "_codon" in stat_axis[stat_index]:
                    plt.ylim((1, 61))
                elif "_aa" in stat_axis[stat_index]:
                    plt.ylim((1, 20))
                elif "at_over_gc" in stat_axis[stat_index]:
                    plt.yscale("log")
                plt.ylabel(stat_names[stat_index], fontsize=24)
                plt.legend(fontsize=12)
            elif 2 == dimension:
                zoom = 3
                import scipy.ndimage
                zommed_x_axis = np.logspace(np.log10(np.min(x_axis)), np.log10(np.max(x_axis)), len(x_axis) * zoom)
                zommed_y_axis = np.logspace(np.log10(np.min(y_axis)), np.log10(np.max(y_axis)), len(y_axis) * zoom)
                x_mesh, y_mesh = np.meshgrid(zommed_x_axis, zommed_y_axis)
                zommed_array = scipy.ndimage.zoom(x_y_array, zoom)
                contourf = plt.contourf(x_mesh, y_mesh, zommed_array, 24, cmap='RdBu',
                                        vmin=np.nanmin(zommed_array),
                                        vmax=np.nanmax(zommed_array))
                contour = plt.contour(contourf,
                                      levels=contourf.levels[::2],
                                      linewidths=(1,), colors='black')
                plt.axis([min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
                cbar = plt.colorbar(contourf)
                cbar.ax.set_ylabel(stat_names[stat_index], fontsize=24)
                # Add the contour line levels to the colorbar
                cbar.add_lines(contour)
                plt.clabel(contour, fmt='%2.1f', colors='black', fontsize=12)
                plt.yscale("log")
                plt.ylabel('${0}$'.format(labels[display_to_order[1]]), fontsize=24)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.xscale("log")
            plt.xlabel('${0}$'.format(labels[display_to_order[0]]), fontsize=24)
            plt.title('${0}={1}$'.format(labels[display_to_order[2]], z_param), fontsize=24)
            plt.tight_layout()
            # plt.savefig("../figures/{0}_{1}_{2}d.svg".format(stat_axis[stat_index], transpose, dimension), format="svg")
            plt.savefig("../figures/{0}_{1}d.png".format(stat_axis[stat_index], dimension), format="png")
            plt.clf()
            plt.close('all')


lambda_axis = np.logspace(np.log10(0.1), np.log10(10), 3)
beta_axis = np.logspace(np.log10(0.1), np.log10(10), 30)
b_axis = [0]
plot_summary_stats(lambda_axis, beta_axis, b_axis, protein="np", dimensions=[1])
