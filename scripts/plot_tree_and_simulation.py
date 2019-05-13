#!python3
import argparse
import numpy as np
import pandas as pd
from ete3 import Tree
from plot_tree import plot_trees_from_traces


def plot_trace(input_simu, input_trace, output_plot):
    simu_dict, color_map_dict = dict(), dict()

    simu_params = {k: v[0] for k, v in pd.read_csv(input_simu + '.parameters.tsv', sep='\t').items()}
    t = Tree(input_simu + ".nhx", format=3)
    root_pop_size = simu_params["population_size"] if "population_size" in simu_params else 1.0
    simu_dict["LogPopulationSize"] = [np.log(float(n.population_size) / root_pop_size) for n in t.traverse()]
    root_age = simu_params["tree_max_distance_to_root_in_year"]

    if "population_size" in simu_params:
        simu_dict["LogMutationRatePerGeneration"] = [np.log(float(n.mutation_rate_per_generation) * root_age) for n in
                                                     t.traverse()]
        simu_dict["LogGenerationTime"] = [np.log(float(n.generation_time)) for n in t.traverse()]
        simu_dict["LogMutationRatePerTime"] = [
            np.log(root_age * float(n.mutation_rate_per_generation) / float(n.generation_time)) for n in t.traverse()]

        simu_dict["LogBranchLength"] = [np.log(n.get_distance(n.up) * float(n.mutation_rate_per_generation) / float(
            n.generation_time)) for n in t.traverse() if not n.is_root()]

        simu_dict["LogTheta"] = [np.log(4 * float(n.mutation_rate_per_generation) * float(n.population_size)) for n in
                                 t.traverse() if n.is_leaf()]
        color_map_dict["LogTheta"] = [np.log(float(n.population_size)) - np.log(float(n.up.population_size)) for n in
                                      t.traverse() if n.is_leaf()]
    else:
        simu_dict["LogMutationRatePerTime"] = [np.log(float(n.mutation_rate) * root_age) for n in t.traverse()]
        simu_dict["LogBranchLength"] = [np.log(n.get_distance(n.up) * float(n.mutation_rate)) for n in t.traverse() if
                                        not n.is_root()]

    simu_dict["BranchTime"] = [n.get_distance(n.up) / root_age for n in t.traverse() if not n.is_root()]

    color_map_branch = [np.log(float(n.population_size)) - np.log(float(n.up.population_size)) for n in t.traverse() if
                        not n.is_root()]
    color_map_dict["BranchTime"] = color_map_branch
    color_map_dict["LogBranchLength"] = color_map_branch

    color_map_nodes = [
        (np.log(float(n.population_size)) - np.log(float(n.up.population_size)) if not n.is_root() else 0.0)
        for n in t.traverse()]
    color_map_dict["LogMutationRatePerGeneration"] = color_map_nodes
    color_map_dict["LogGenerationTime"] = color_map_nodes
    color_map_dict["LogMutationRatePerTime"] = color_map_nodes

    plot_trees_from_traces(input_trace, output_plot, simu_dict, color_map_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--simulation', required=True, type=str, dest="simulation")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, nargs='+', dest="trace")
    args = parser.parse_args()
    plot_trace(args.simulation, args.trace, args.output)
