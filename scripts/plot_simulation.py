#!python3
from ete3 import Tree
from plot_module import plot_correlation, plot_tree, convertible_to_float
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import os

min_max_annot = ()


def plot_simulation(input_simu, args_output):
    t = Tree(input_simu, format=1)

    args_nodes = set()
    for node in t.traverse():
        args_nodes = args_nodes.union(node.features)

    branch_dict = {}
    for arg in args_nodes:
        if arg != "dist" and arg != "support":
            values = np.array([float(getattr(n, arg)) for n in t.traverse() if
                               arg in n.features and convertible_to_float(getattr(n, arg))])
            if len(values) > 1 and len(values) == len(list(t.traverse())):
                root_pop_size = float(getattr(t.get_tree_root(), arg))
                for n in t.traverse():
                    n.add_feature("Log" + arg, np.log(float(getattr(n, arg)) / root_pop_size))
                plot_tree(t, "Log" + arg, "{0}/tree.{1}.png".format(args_output, arg), min_max_annot=min_max_annot)
            if len(values) > 1 and ("Branch" in arg) and (("dNd" in arg) or ("LogNe" in arg)):
                branch_dict[arg] = values

    der_pop_size = [(np.log10(float(n.population_size)) - float(n.Branch_LogNe)) for n in t.traverse() if
                    not n.is_root()]
    plot_correlation("{0}/correlation.Ne.dNdS.svg".format(args_output), branch_dict, {}, der_pop_size,
                     min_max_annot=[0, 1.0], global_xy=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, default='',
                        dest="t", metavar="<tree>", help="The tree to be drawn")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    args = parser.parse_args()
    plot_simulation(args.t, args.output)
