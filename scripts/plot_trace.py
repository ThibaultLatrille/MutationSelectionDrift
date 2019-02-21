#!python3
import argparse
import os
from csv import reader
from ete3 import Tree, TreeStyle, TextFace, faces, CircleFace
import numpy as np
import statsmodels.api as sm
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

my_dpi = 128
RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"


def tex_f(f):
    if 1 <= f < 10:
        return "{0:.2g}".format(f)
    else:
        return "{0:.2e}".format(f)


def to_float(element):
    try:
        return float(element)
    except ValueError:
        return element


def open_tsv(filepath):
    if os.path.exists(filepath):
        with open(filepath, 'r') as tsv_open:
            tsv = list(reader(tsv_open, delimiter='\t'))
        return tsv
    else:
        return [[]]


def is_ultrametric(tree):
    d = [tree.get_distance(leaf) for leaf in tree.iter_leaves()]
    if abs(max(d) - min(d)) / max(d) < 1e-4:
        print("The tree is ultrametric.")
        return True
    else:
        print("The tree is not ultrametric.")
        print("Max distance from root to leaves: {0}.".format(max(d)))
        print("Min distance from root to leaves: {0}.".format(min(d)))
        return False


def layout(node, arg, min_arg, max_arg, filename, columns):
    if arg + "." + filename in node.features:
        if min_arg == max_arg:
            radius = 15
        else:
            radius = 15 * (float(getattr(node, arg + "." + filename)) - min_arg) / (max_arg - min_arg) + 5
        circle = CircleFace(radius=radius, color="RoyalBlue", style="sphere")
        circle.opacity = 0.3
        faces.add_face_to_node(circle, node, 0, position="float")
        for col, attr in enumerate(columns):
            faces.add_face_to_node(TextFace(tex_f(float(getattr(node, arg + "." + attr))) + " "), node, col,
                                   position="aligned")


def plot_trace(input_simu, input_trace, output_plot, burn_in):
    simu_params, traces = dict(), dict()

    params = open_tsv(input_simu + '.parameters.tsv')
    for i, param in enumerate(params[0]):
        simu_params[param] = to_float(params[1][i])

    filenames = ["Simulation", "Watterson_Simulation", "WattersonSynonymous_Simulation", "Pairwise_Simulation",
                 "PairwiseSynonymous_Simulation"]
    for filepath in input_trace:
        filenames.append(os.path.basename(filepath))
        trace = open_tsv(filepath + '.trace')
        for i, param in enumerate(trace[0]):
            if param not in traces:
                traces[param] = dict()
            traces[param][os.path.basename(filepath)] = [float(line[i]) for line in trace[(1 + burn_in):]]

    tree = Tree(input_simu + ".nhx", format=3)
    is_ultrametric(tree)
    for node in tree.traverse():
        rate = max(float(node.mutation_rate_per_generation) / float(node.generation_time_in_year), 1e-30)
        node.add_feature("LogNe.Simulation", np.log(float(node.population_size)))
        node.add_feature("LogMutRate.Simulation", np.log(rate))
        if not node.is_root():
            time = float(node.get_distance(node.up)) / simu_params["tree_max_distance_to_root_in_year"]
            node.add_feature("Time.Simulation", float(time))
            node.add_feature("LogLength.Simulation", np.log(time * rate))
            node.add_feature("dNdS.Simulation", float(node.dNdS))
        if node.is_leaf():
            node.add_feature("Theta.Simulation", 4 * float(node.population_size) * float(node.mutation_rate_per_generation))
            node.add_feature("Theta.Watterson_Simulation", float(node.theta_watterson))
            node.add_feature("Theta.WattersonSynonymous_Simulation", float(node.theta_watterson_syn))
            node.add_feature("Theta.Pairwise_Simulation", float(node.theta_pairwise))
            node.add_feature("Theta.PairwiseSynonymous_Simulation", float(node.theta_pairwise_syn))
    for param, traces_param in traces.items():
        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        for filename, param_trace in sorted(traces_param.items(), key=lambda x: x[0]):
            style = "-"
            if "False" in filename:
                style = "--"
            plt.plot(range(len(param_trace)), param_trace, style, alpha=0.5, linewidth=1, label=filename)

            if param[0] == "*":
                node_name = param.split("_")[-1]

                if node_name != "":
                    node = next(tree.iter_search_nodes(name=node_name))
                else:
                    node = tree.get_tree_root()

                if "*NodePopSize" in param:
                    node.add_feature("LogNe." + filename, np.mean(param_trace))
                    node.add_feature("LogNe_var." + filename, np.var(param_trace))
                elif "*NodeRate" in param:
                    node.add_feature("LogMutRate." + filename, np.mean(param_trace))
                    node.add_feature("LogMutRate_var." + filename, np.var(param_trace))
                elif "*BranchTime" in param:
                    assert (not node.is_root())
                    node.add_feature("Time." + filename, np.mean(param_trace))
                    node.add_feature("Time_var." + filename, np.var(param_trace))
                elif "*BranchLength" in param:
                    assert (not node.is_root())
                    node.add_feature("LogLength." + filename, np.mean(np.log(param_trace)))
                    node.add_feature("LogLength_var." + filename, np.var(np.log(param_trace)))
                elif "*BranchdNdS" in param:
                    assert (not node.is_root())
                    node.add_feature("dNdS." + filename, np.mean(param_trace))
                    node.add_feature("dNdS_var." + filename, np.var(param_trace))
                elif "*Theta" in param:
                    assert (not node.is_leaf())
                    node.add_feature("Theta." + filename, np.mean(param_trace))
                    node.add_feature("Theta_var." + filename, np.var(param_trace))
        if param.lower() in [s.lower() for s in simu_params]:
            # If it can be found in the Simulation parameters
            plt.axhline(y=simu_params[param], xmin=0.0, xmax=1.0, color='r', label="Simulation")

        plt.xlabel('Point')
        plt.ylabel(param)
        plt.legend()
        plt.tight_layout()
        plt.savefig('{0}/trace.{1}.png'.format(output_plot, param), format='png')
        plt.clf()
        plt.close('all')

    for arg in ["Time", "Theta", "dNdS", "LogLength", "LogMutRate", "LogNe"]:
        axis_dict = dict()
        err_dict = dict()
        for filename in filenames:
            attr = arg + "." + filename
            var = arg + "_var." + filename
            axis = np.array([getattr(node, attr) for node in tree.traverse() if attr in node.features])
            err = np.array([getattr(node, var) for node in tree.traverse() if var in node.features])
            if len(axis) > 0:
                axis_dict[filename] = axis
                if len(err) > 0:
                    err_dict[filename] = np.sqrt(err)
                else:
                    err_dict[filename] = np.array([0] * len(axis))
        for filename, axis in axis_dict.items():
            max_arg, min_arg = max(axis), min(axis)
            ts = TreeStyle()
            ts.show_leaf_name = True
            ts.complete_branch_lines_when_necessary = False
            columns = sorted(axis_dict)
            ts.layout_fn = lambda x: layout(x, arg, min_arg, max_arg, filename, columns)
            for col, name in enumerate(columns):
                nameF = TextFace(name, fsize=7)
                nameF.rotation = -90
                ts.aligned_header.add_face(nameF, column=col)
            ts.title.add_face(TextFace("{1} in {2}".format(output_plot, arg, filename), fsize=20), column=0)
            tree.render("{0}/nhx.{1}.{2}.png".format(output_plot, arg, filename), tree_style=ts)

        f, axs = plt.subplots(len(axis_dict), len(axis_dict),
                              figsize=(len(axis_dict) * 640 / my_dpi, len(axis_dict) * 480 / my_dpi), dpi=my_dpi)

        if len(axis_dict) == 1:
            axs = [[axs]]

        def min_max(axis, err):
            eps = 0.05
            min_axis, max_axis = min(axis - err), max(axis + err)
            min_axis -= (max_axis - min_axis) * eps
            max_axis += (max_axis - min_axis) * eps
            if min_axis == max_axis:
                return min_axis - eps, max_axis + eps
            else:
                return min_axis, max_axis

        for row, (row_filename, row_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
            for col, (col_filename, col_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
                ax = axs[row][col]
                ax.errorbar(col_axis, row_axis,
                            xerr=err_dict[col_filename], yerr=err_dict[row_filename],
                            fmt='o', color=BLUE, ecolor=GREEN, label=r"${0}$ nodes".format(len(row_axis)))

                min_col, max_col = min_max(col_axis, err_dict[col_filename])
                idf = np.linspace(min_col, max_col, 30)
                if col == 0:
                    ax.set_ylabel(row_filename)
                if row == len(axis_dict) - 1:
                    ax.set_xlabel(col_filename)
                ax.set_xlim((min_col, max_col))
                min_row, max_row = min_max(row_axis, err_dict[row_filename])
                ax.set_ylim((min_row, max_row))
                if row_filename != col_filename and len(set(col_axis)) > 1 and len(set(row_axis)) > 1:
                    model = sm.OLS(row_axis, sm.add_constant(col_axis))
                    results = model.fit()
                    b, a = results.params[0:2]
                    ax.plot(idf, a * idf + b, '-', color=RED, label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
                        float(a), abs(float(b)), results.rsquared, "+" if float(b) > 0 else "-"))
                    ax.legend()
        plt.tight_layout()
        plt.savefig('{0}/correlation.{1}.png'.format(output_plot, arg), format='png')
        plt.clf()
        plt.close('all')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--simu', required=True, type=str, dest="simu")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, nargs='+', dest="trace")
    parser.add_argument('-b', '--burn_in', required=False, type=int, default=0, dest="burn_in")
    args = parser.parse_args()
    plot_trace(args.simu, args.trace, args.output, args.burn_in)
