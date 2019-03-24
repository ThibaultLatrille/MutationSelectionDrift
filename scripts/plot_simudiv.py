#!python3
from plot_module import *
import argparse


def plot_simulation(input_simu):
    t = Tree(input_simu, format=3)

    args_nodes = set()
    for node in t.traverse():
        args_nodes = args_nodes.union(node.features)

    branch_dict = {}
    for arg in args_nodes:
        values = [float(getattr(n, arg)) for n in t.traverse() if
                  arg in n.features and convertible_to_float(getattr(n, arg))]
        if len(values) > 1:
            if ("Branch" in arg) and (("dNd" in arg) or ("LogNe" in arg)):
                branch_dict[arg] = values

    plot_correlation("{0}.correlation.png".format(input_simu), branch_dict, {})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, default='',
                        dest="t", metavar="<tree>", help="The tree to be drawn")
    args = parser.parse_args()
    plot_simulation(args.t)