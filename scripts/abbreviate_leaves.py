#!python3
from ete3 import Tree
import argparse


def tree_plot(input_tree):
    t = Tree(input_tree, format=1)
    names = set()
    for leaf in t:
        name = leaf.name[0].upper() + leaf.name.split("_")[1][:3].lower()
        assert(name not in names)
        names.add(name)
        print(name)
        leaf.name = name
    t.write(format=1, outfile=input_tree + ".abbr")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str,
                        default='../DataEmpirical/sp73_OrthoMam/rootedtree.nhx', dest="t", metavar="<tree>",
                        help="The tree to be re-written")
    args = parser.parse_args()
    tree_plot(args.t)
