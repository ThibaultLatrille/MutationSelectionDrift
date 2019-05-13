#!python3
import os
import argparse
import numpy as np
from glob import glob
from ete3 import Tree
from plot_module import plot_correlation, plot_tree


def plot_trees_from_traces(input_trace, output_plot, simu_dict, color_map_dict):
    axis_trees, axis_filenames = dict(), dict()

    for filepath in input_trace:
        for tree_path in sorted(glob("{0}.*.nhx".format(filepath))):
            feature = tree_path.replace(filepath + ".", "").replace(".nhx", "")
            filename = os.path.basename(filepath)
            tree = Tree(tree_path, format=1)

            if feature not in axis_trees:
                axis_trees[feature] = []
            if feature not in axis_filenames:
                axis_filenames[feature] = []

            axis_filenames[feature].append(filename)
            axis_trees[feature].append(tree)
            if len([n for n in tree.traverse() if feature in n.features]) == len(list(tree.traverse())):
                plot_tree(tree.copy(), feature, "{0}/{1}.{2}.svg".format(output_plot, filename, feature))

    for feature in axis_trees:
        axis_dict, err_dict = dict(), dict()
        if feature in simu_dict:
            axis_dict["Simulation"] = simu_dict[feature]
        for filename, tree in zip(axis_filenames[feature], axis_trees[feature]):
            values = np.array([float(getattr(n, feature)) for n in tree.traverse() if feature in n.features])
            min_values = np.array(
                [float(getattr(n, feature + "_min")) for n in tree.traverse() if feature + "_min" in n.features])
            max_values = np.array(
                [float(getattr(n, feature + "_max")) for n in tree.traverse() if feature + "_max" in n.features])
            axis_dict[filename] = values
            err_dict[filename] = np.vstack((np.abs(values - min_values), np.abs(max_values - values)))

        if len(axis_dict) > 1:
            path = '{0}/correlation.{1}.png'.format(output_plot, feature)
            if feature in color_map_dict:
                plot_correlation(path, axis_dict, err_dict, color_map_dict[feature])
            else:
                plot_correlation(path, axis_dict, err_dict, [])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, nargs='+', dest="trace")
    args = parser.parse_args()
    plot_trees_from_traces(args.trace, args.output, {}, {})
