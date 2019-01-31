# GLOBAL IMPORTS
from analysis import *
import os
import matplotlib.pyplot as plt
from matplotlib import markers
markers_array = list(markers.MarkerStyle.markers.keys())

protein = "np"
inference = "0-3-1"
title = "Inference using mean-field parametric omega"
current_dir = "/home/thibault/SimuEvol/simulated"
for file in sorted(os.listdir(current_dir + "/data_prefs")):
    prefix = file.strip().replace(".txt", "")
    hyphy_path = "{0}/data_hyphy/{1}".format(current_dir, prefix)
    if os.path.isdir(hyphy_path) and (protein in prefix):
        nested_dict = nested_dict_init()
        protein_prefs_path = "{0}/data_prefs/{1}".format(current_dir, file)
        batchfiles = [batch[:-3] for batch in os.listdir(hyphy_path) if ((batch[-3:] == ".bf") and (inference in batch))]

        for batch in batchfiles:
            at_gc_pct = float(batch.split("_m")[-1].split("_")[0])
            simu_evol_path = "{0}/data_alignment/{1}_{2}".format(current_dir, prefix, "_".join(batch.split("_")[:-1]))
            simuevol_dico = dico_from_file(simu_evol_path + ".txt")

            nuc_freqs, nbr_sites, nbr_species = extract_nuc_pct(simu_evol_path + ".fasta")
            at_gc_pct_obs = (nuc_freqs['A'] + nuc_freqs['T']) / (nuc_freqs['G'] + nuc_freqs['C'])
            nested_dict["AT/GC_obs"][at_gc_pct] = at_gc_pct_obs
            nested_dict["w_obs"][at_gc_pct] = simuevol_dico["w"]

            hyphy_result = "{0}/{1}.bf_hyout.txt".format(hyphy_path, batch)
            hyphy_dico = dico_from_file(hyphy_result)
            hyphy_dico["n"] = nbr_sites
            format_hyphy_dico(hyphy_dico)

            nested_dict["w_inf"][at_gc_pct] = hyphy_dico["w"]

            if ("pnG" in hyphy_dico) and ("pnG" in hyphy_dico):
                gc_pct = hyphy_dico["pnG"] + hyphy_dico["pnC"]
                nested_dict["lambda_inf"][at_gc_pct] = (1 - gc_pct) / gc_pct

            nested_dict["AT/GC_inf"][at_gc_pct] = equilibrium_lambda(hyphy_dico)

        my_dpi = 96
        fig, ax = plt.subplots()
        index = 0
        list_plot = list()
        list_plot.append({"experiment": "AT/GC_obs", "color": BLUE, "linestyle": '-', "linewidth": 2, "label": "$AT/GC$ observed"})
        list_plot.append({"experiment": "AT/GC_inf", "color": YELLOW, "linestyle": '--', "linewidth": 4, "label": "$\widehat{AT/GC}$ inferred"})
        list_plot.append({"experiment": "lambda_inf", "color": GREEN, "linestyle": '--', "linewidth": 4, "label": "$\widehat{\lambda}$ inferred"})

        x_list = sorted(nested_dict["AT/GC_obs"].keys())
        ax.plot(x_list, x_list, color="black", linestyle='-', linewidth=2, label="y=x")

        omega_obs = np.array([nested_dict["w_obs"][k] for k in x_list])
        omega_inf = np.array([nested_dict["w_inf"][k] for k in x_list])
        print("|w_obs-w_inf|/w_obs = {0:.2f}%".format(100 * np.mean(np.abs((omega_inf - omega_obs)) / omega_obs)))

        for param in list_plot:
            x_list = sorted(nested_dict[param["experiment"]].keys())
            y_list = [nested_dict[param["experiment"]][k] for k in x_list]
            ax.plot(x_list, y_list, linestyle=param["linestyle"], label=param["label"],
                    color=param["color"], linewidth=param["linewidth"])

        lambda_obs = np.array(x_list)
        lambda_inf = np.array([nested_dict["lambda_inf"][k] for k in x_list])
        print("|lambda_obs-lambda_inf|/lambda_obs = {0:.2f}%".format(100 * np.mean(np.abs((lambda_inf - lambda_obs)) / lambda_obs)))

        ax.set_xscale('log')
        ax.set_xlabel('$\lambda$ used for the simulation')
        ax.set_yscale('log')
        ax.set_ylabel('$\lambda$ inferred')
        ax.legend()
        ax.set_title(title)
        plt.tight_layout()
        plt.savefig("../figures/np_simulated_{0}.svg".format(inference), format="svg")
        plt.savefig("../figures/np_simulated_{0}.png".format(inference), format="png")
        plt.show()
        plt.clf()
        plt.close('all')
