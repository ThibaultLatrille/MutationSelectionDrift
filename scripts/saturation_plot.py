# GLOBAL IMPORTS
from analysis import *
import matplotlib.pyplot as plt
import os

current_dir = "/home/thibault/SimuEvol/"
nbr_steps = 14
nbr_points = 5
mut_bias_list = np.logspace(-1, 1, nbr_points)
x_axis = 1000 * np.logspace(-1, 1, nbr_steps)
y_axis = np.logspace(-1.5, 1.5, nbr_steps)
y, x = np.meshgrid(y_axis, x_axis)
for mut_bias_id in range(nbr_points):
    omega_sim = np.zeros((len(x), len(y)))
    omega_mg = np.zeros((len(x), len(y)))
    omega_proj = np.zeros((len(x), len(y)))
    omega_proj_eps = np.zeros((len(x), len(y)))

    for x_i, nbr_sites in enumerate(x_axis):
        for y_i, mu in enumerate(y_axis):
            protein = "np"
            mixture = True
            alpha = 0.1
            nbr_sites = int(nbr_sites)
            prefix = "{0}_{1}_{2}_{3}_{4}".format(nbr_sites, protein, mixture, alpha, mu)

            hyphy_path = "{0}/data_hyphy/{1}".format(current_dir, prefix)
            if os.path.isdir(hyphy_path):
                protein_prefs_path = "{0}/data_prefs/{1}.txt".format(current_dir, prefix)

                for qsub_id, mut_bias in enumerate(mut_bias_list):
                    if qsub_id == mut_bias_id:
                        id_sufix = "id{0}_m{1:.5f}".format(qsub_id, mut_bias)
                        name = "{0}_{1}".format(prefix, id_sufix)
                        ali_path = "{0}/data_alignment/{1}".format(current_dir, name)

                        simuevol_dico = dico_from_file(ali_path + ".txt")
                        if abs(simuevol_dico["w"] - simuevol_dico["w0"]) <= 0.2 and (simuevol_dico["w"] != 0.0):
                            omega_sim[x_i, y_i] = simuevol_dico["w"]
                        else:
                            omega_sim[x_i, y_i] = None

                        rate = 0  # can be 0, 1, or 5
                        freq = 3  # can be 1 or 3
                        for omega in [1, 20, 95]:
                            param = "{0}-{1}-{2}".format(rate, freq, omega)
                            hyphy_result = "{0}/{1}_{2}.bf_hyout.txt".format(hyphy_path, id_sufix, param)
                            if omega == 1:
                                matrix = omega_mg
                            elif omega == 20:
                                matrix = omega_proj_eps
                            else:
                                matrix = omega_proj

                            if os.path.isfile(hyphy_result) and (omega_sim[x_i, y_i] is not None):
                                hyphy_dico = dico_from_file(hyphy_result)
                                hyphy_dico["n"] = nbr_sites
                                format_hyphy_dico(hyphy_dico)
                                matrix[x_i, y_i] = omega_sim[x_i, y_i] / hyphy_dico["w"]
                            else:
                                matrix[x_i, y_i] = None

    titles = ['$\omega_{sim}$',
              '$\omega_{sim} / \hat{\omega}_{MG}$',
              '$\omega_{sim} / \hat{\omega}_{95}$',
              '$\omega_{sim} / \hat{\omega}_{20}$']
    m_list = [omega_mg, omega_proj, omega_proj_eps]
    mini = min([np.nanmin(m) for m in m_list])
    maxi = min([np.nanmax(m) for m in m_list])
    for index, z in enumerate([omega_sim] + m_list):
        plt.subplot(2, 2, index + 1)
        if index == 0:
            vmin, vmax = np.nanmin(z), np.nanmax(z)
        else:
            vmin, vmax = mini, maxi
        plt.pcolor(x, y, z, cmap='RdBu', vmin=vmin, vmax=vmax)
        plt.xscale("log")
        plt.yscale("log")
        plt.title(titles[index])
        plt.xlabel('Nbr of sites')
        plt.ylabel('$\mu$')
        # set the limits of the plot to the limits of the data
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        plt.colorbar()
    plt.suptitle('$\lambda = {0:.2f}$'.format(mut_bias_list[mut_bias_id]))
    plt.show()
