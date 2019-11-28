#!python3
import os
import argparse
import numpy as np
import pandas as pd
from glob import glob
from ete3 import Tree
from plot_module import plot_correlation, plot_tree, to_float

min_max_annot = ()


def remove_units(str):
    return str.replace("(g)", "").replace("(days)", "").replace("(yrs)", "").replace("(kg)", "").replace("(cm)", "")


def plot_trees_from_traces(input_trace, output_plot, simu_dict, color_map_dict, simu_tree):
    axis_trees, axis_filenames = dict(), dict()

    for filepath in input_trace:
        for tree_path in sorted(glob("{0}.*.nhx".format(filepath))):
            feature = remove_units(tree_path.replace(filepath + ".", "").replace(".nhx", ""))
            filename = os.path.basename(filepath)
            with open(tree_path, 'r') as tree_file:
                tree_str = remove_units(tree_file.readline())
                if tree_str.count("-nan") > 0:
                    continue
                tree = Tree(tree_str, format=1)

            if simu_tree:
                for n_inf, n_simu in zip(tree.traverse(), simu_tree.traverse()):
                    assert (sorted(n_simu.get_leaf_names()) == sorted(n_inf.get_leaf_names()))

            if feature not in axis_trees:
                axis_trees[feature] = []
            if feature not in axis_filenames:
                axis_filenames[feature] = []

            axis_filenames[feature].append(filename)
            axis_trees[feature].append(tree)
            if len([n for n in tree.traverse() if feature in n.features]) == len(list(tree.traverse())):
                plot_tree(tree.copy(), feature, "{0}/{1}.{2}.png".format(output_plot, filename, feature), min_max_annot=min_max_annot)

    for feature in axis_trees:
        axis_dict, err_dict = dict(), dict()
        if feature in simu_dict:
            axis_dict["Simulation"] = simu_dict[feature]
        for filename, tree in zip(axis_filenames[feature], axis_trees[feature]):
            values = np.array([to_float(getattr(n, feature)) for n in tree.traverse() if feature in n.features])
            min_values = np.array(
                [to_float(getattr(n, feature + "_min")) for n in tree.traverse() if feature + "_min" in n.features])
            max_values = np.array(
                [to_float(getattr(n, feature + "_max")) for n in tree.traverse() if feature + "_max" in n.features])
            axis_dict[filename] = values
            err_dict[filename] = np.vstack((np.abs(values - min_values), np.abs(max_values - values)))

        if len(axis_dict) > 1:
            path = '{0}/correlation.{1}.svg'.format(output_plot, feature)

            if feature in color_map_dict:
                plot_correlation(path, axis_dict, err_dict, color_map_dict[feature], min_max_annot=min_max_annot,
                                 global_xy=False)
            else:
                plot_correlation(path, axis_dict, err_dict, [], min_max_annot=(), global_xy=False)


def open_simulation(input_simu):
    simu_dict, color_map_dict = dict(), dict()

    simu_params = {k: v[0] for k, v in pd.read_csv(input_simu + '.parameters.tsv', sep='\t').items()}
    t = Tree(input_simu + ".nhx", format=1)
    root_pop_size = float(t.population_size)
    simu_dict["LogPopulationSize"] = [np.log(float(n.population_size) / root_pop_size) for n in t.traverse()]
    root_age = simu_params["tree_max_distance_to_root_in_year"]

    simu_dict["ContrastPopulationSize"] = [
        (np.log(float(n.population_size)) - np.log(float(n.up.population_size))) / np.sqrt(
            n.get_distance(n.up) / root_age) for n in t.traverse() if not n.is_root()]

    if "population_size" in simu_params:
        simu_dict["LogMutationRatePerGeneration"] = [np.log(float(n.mutation_rate_per_generation) * root_age) for n in
                                                     t.traverse()]
        simu_dict["LogGenerationTime"] = [np.log(float(n.generation_time)) for n in t.traverse()]
        simu_dict["LogMutationRatePerTime"] = [
            np.log(root_age * float(n.mutation_rate_per_generation) / float(n.generation_time)) for n in t.traverse()]

        simu_dict["Log10Theta"] = [np.log10(4 * float(n.mutation_rate_per_generation) * float(n.population_size)) for n
                                   in
                                   t.traverse() if n.is_leaf()]
        color_map_dict["Log10Theta"] = [np.log(float(n.population_size) / root_pop_size) for n in t.traverse() if
                                        n.is_leaf()]
    else:
        simu_dict["LogMutationRatePerTime"] = [np.log(float(n.mutation_rate) * root_age) for n in t.traverse()]

    simu_dict["BranchTime"] = [n.get_distance(n.up) / root_age for n in t.traverse() if not n.is_root()]
    simu_dict["Log10BranchLength"] = [
        np.log10(n.get_distance(n.up) * float(n.Branch_mutation_rate_per_generation) / float(n.Branch_generation_time))
        for n in t.traverse() if not n.is_root()]

    color_map_branch = [np.log(float(n.population_size) / root_pop_size) for n in t.traverse() if not n.is_root()]
    color_map_dict["BranchTime"] = color_map_branch
    color_map_dict["Log10BranchLength"] = color_map_branch
    color_map_dict["ContrastPopulationSize"] = color_map_branch

    color_map_nodes = [np.log(float(n.population_size) / root_pop_size) for n in t.traverse()]
    color_map_dict["LogMutationRatePerGeneration"] = color_map_nodes
    color_map_dict["LogGenerationTime"] = color_map_nodes
    color_map_dict["LogMutationRatePerTime"] = color_map_nodes
    color_map_dict["LogPopulationSize"] = color_map_nodes
    return simu_dict, color_map_dict, t


def open_tsv_population_size(tree_file, tsv_file):
    t = Tree(tree_file, format=1)
    csv = pd.read_csv(tsv_file, header=None, sep='\t')
    for index, (leaf_1, leaf_2, _, ne, _) in csv.iterrows():
        if leaf_1 == leaf_2:
            leaves = t.get_leaves_by_name(leaf_1)
            assert (len(leaves) == 1)
            n = leaves[0]
        else:
            n = t.get_common_ancestor([leaf_1, leaf_2])
        n.pop_size = ne

    pop_size_dict = dict()
    root_pop_size = float(t.pop_size)
    pop_size_dict["LogPopulationSize"] = [np.log(float(n.pop_size) / root_pop_size) for n in t.traverse()]
    return pop_size_dict, t


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--simulation', required=False, default="", type=str, dest="simulation")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, nargs='+', dest="trace")
    parser.add_argument('--tree', required=False, default="", type=str, dest="tree")
    parser.add_argument('--tsv', required=False, default="", type=str, dest="tsv")

    args = parser.parse_args()
    if args.simulation != "":
        simu, color_map, simu_tree = open_simulation(args.simulation)
        plot_trees_from_traces(args.trace, args.output, simu, color_map, simu_tree)
    elif args.tree != "" and args.tsv != "":
        pop_size, input_tree = open_tsv_population_size(args.tree, args.tsv)
        plot_trees_from_traces(args.trace, args.output, pop_size, {}, input_tree)
    else:
        plot_trees_from_traces(args.trace, args.output, {}, {}, False)
