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
plot = False


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
    simu_params, simu_nodes, simu_taxa, traces = dict(), dict(), dict(), dict()

    params = open_tsv(input_simu + '.parameters.tsv')
    for i, param in enumerate(params[0]):
        simu_params[param] = to_float(params[1][i])

    nodes = open_tsv(input_simu + '.nodes.tsv')
    header = nodes[0]
    index = {v: k for k, v in enumerate(header)}["index"]
    for line in nodes[1:]:
        simu_nodes[line[index]] = dict()
        for i, value in enumerate(line):
            simu_nodes[line[index]][header[i]] = to_float(value)

    taxa = open_tsv(input_simu + '.tsv')
    header = taxa[0]
    index = {v: k for k, v in enumerate(header)}["taxon_name"]
    for line in nodes[1:]:
        simu_nodes[line[index]] = dict()
        for i, value in enumerate(line):
            simu_nodes[line[index]][header[i]] = to_float(value)

    filenames = ["simulation"]
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
        if not node.is_root():
            time = float(node.get_distance(node.up)) / simu_params["tree_max_distance_to_root_in_year"]
            node.add_feature("Time.simulation", time)
            node.add_feature("Length.simulation", np.log(time * rate))
        node.add_feature("LogNe.simulation", np.log(float(node.population_size)))
        node.add_feature("LogMutRate.simulation", np.log(rate))
        if node.is_leaf():
            node.add_feature("LogTheta.simulation",
                             np.log(4 * float(node.population_size) * float(node.mutation_rate_per_generation)))

    for param, traces_param in traces.items():
        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        for filename, param_trace in sorted(traces_param.items(), key=lambda x: x[0]):
            style = "-"
            if "False" in filename:
                style = "--"
            if ("NodePopSize" in param) and (max(param_trace) > 10):
                print("Issue with " + filename + ", " + param)
                y_values = [min(i, 10) for i in param_trace]
                plt.plot(range(len(param_trace)), y_values, style, alpha=0.5, linewidth=1, label="!" + filename)
            else:
                plt.plot(range(len(param_trace)), param_trace, style, alpha=0.5, linewidth=1, label=filename)

            if param[0] == "*" or param[0] == "@":
                node_name = param.split("_")[-1]

                if node_name != "":
                    node = next(tree.iter_search_nodes(name=node_name))
                else:
                    node = tree.get_tree_root()

                if "*NodePopSize" in param:
                    node.add_feature("LogNe." + filename, np.mean(param_trace))
                elif "*NodeRate" in param:
                    node.add_feature("LogMutRate." + filename, np.mean(param_trace))
                elif "*BranchTime" in param:
                    assert (not node.is_root())
                    node.add_feature("Time." + filename, np.mean(param_trace))
                elif "*BranchLength" in param:
                    assert (not node.is_root())
                    node.add_feature("Length." + filename, np.log(np.mean(param_trace)))
                elif "*Theta" in param:
                    node.add_feature("LogTheta." + filename, np.log(np.mean(param_trace)))
        if param.lower() in [s.lower() for s in simu_params]:
            # If it can be found in the simulation parameters
            plt.axhline(y=simu_params[param], xmin=0.0, xmax=1.0, color='r', label="Simulation")

        plt.xlabel('Point')
        plt.ylabel(param)
        plt.legend()
        plt.tight_layout()
        plt.savefig('{0}/trace.{1}.png'.format(output_plot, param), format='png')
        plt.clf()
        plt.close('all')

    for arg in ["Time", "Length", "LogTheta", "LogMutRate", "LogNe"]:
        axis_dict = dict()
        for filename in filenames:
            attr = arg + "." + filename
            axis = [float(getattr(node, attr)) for node in tree.traverse() if attr in node.features]
            if len(axis) > 0:
                axis_dict[filename] = axis

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

        def min_max(axis):
            eps = 0.05
            min_axis, max_axis = min(axis), max(axis)
            min_axis -= (max_axis - min_axis) * eps
            max_axis += (max_axis - min_axis) * eps
            if min_axis == max_axis:
                return min_axis - eps, max_axis + eps
            else:
                return min_axis, max_axis

        for row, (row_filename, row_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
            for col, (col_filename, col_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
                ax = axs[row][col]
                ax.scatter(col_axis, row_axis, linewidth=3, label=r"${0}$ nodes".format(len(row_axis)))

                min_col, max_col = min_max(col_axis)
                idf = np.linspace(min_col, max_col, 30)
                if col == 0:
                    ax.set_ylabel(row_filename)
                if row == len(axis_dict) - 1:
                    ax.set_xlabel(col_filename)
                ax.set_xlim((min_col, max_col))
                min_row, max_row = min_max(row_axis)
                ax.set_ylim((min_row, max_row))
                if row_filename != col_filename and len(set(col_axis)) > 1 and len(set(row_axis)) > 1:
                    model = sm.OLS(row_axis, sm.add_constant(col_axis))
                    results = model.fit()
                    b, a = results.params[0:2]
                    ax.plot(idf, a * idf + b, 'r-', label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
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
