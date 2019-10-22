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

    nodes_file = input_simu.replace(".nhx", ".nodes.tsv")
    if os.path.isfile(nodes_file):
        nodes = pd.read_csv(nodes_file, sep='\t')
        ss_nodes = ["Branch_dN", "Branch_dS", "Branch_dNdS"]
        leaves = pd.read_csv(input_simu.replace(".nhx", ".leaves.tsv"), sep='\t')
        ss_leaves = []
        for method in ["pairwise", "watterson"]:
            for subset in ["", "_syn", "_non_syn"]:
                ss_leaves.append("theta_{0}{1}".format(method, subset))
        data = {"nodes": (nodes, ss_nodes), "leaves": (leaves, ss_leaves)}

        for csv, ss_list in data.values():
            for ss in ss_list:
                cols, groups = zip(*csv.groupby(by="taxon_name"))
                values = np.array([g[ss].values for g in groups]).T
                df = pd.DataFrame(values, columns=cols)
                df.boxplot(column=list(cols))
                for key in [ss, ss + "_pred", ss + "_flow"]:
                    values = [float(getattr(n, key)) for n in t.traverse() if key in n.features]
                    if len(values) > 0 and len(values) == len(cols):
                        plt.scatter(list(cols), values, label="Exome " + key)

                if "theta" in ss:
                    plt.yscale("log")
                plt.ylabel(ss)
                plt.xticks(rotation='vertical')
                plt.legend()
                plt.savefig("{0}/boxplot.{1}.png".format(args_output, ss), bbox_inches='tight')
                plt.clf()
                plt.close("all")

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
