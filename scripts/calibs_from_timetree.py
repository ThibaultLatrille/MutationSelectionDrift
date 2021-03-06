#!python3
import os
import argparse
import pandas as pd
from ete3 import Tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--nwk', default="../DataEmpirical/PrimatesBinaryLHTShort/rootedtree.nwk",
                        required=False, type=str, dest="nwk")
    parser.add_argument('-t', '--tree', default="../DataEmpirical/PrimatesBinaryLHTShort/rootedtree.nwk.annotated",
                        required=False, type=str, dest="tree")
    args = parser.parse_args()
    tree = Tree(args.tree, format=1)
    nwk = Tree(args.nwk, format=1)
    assert (len(tree) == len(nwk))
    tree_leaf_names = set(tree.get_leaf_names())
    for leaf in nwk:
        if leaf.name not in tree_leaf_names:
            try:
                leaf.name = next(n for n in tree_leaf_names if leaf.name.split("_")[0] == n.split("_")[0])
            except StopIteration:
                leaf.name = next(n for n in tree_leaf_names if leaf.name.split("_")[1] == n.split("_")[1])
    assert(set(nwk.get_leaf_names()) == tree_leaf_names)
    root_age = nwk.get_closest_leaf()[1]
    df = [["Root", root_age, root_age, root_age]]
    for n in nwk.iter_descendants(strategy='postorder'):
        if not n.is_leaf():
            name = tree.get_common_ancestor(n.get_leaf_names()).name
            if n.dist <= 0:
                print("Node " + name + " is attached to it's parent, thus it's discarded.")
                continue
            age = n.get_closest_leaf()[1]
            df += [[name, age, age, age]]

    header = ["NodeName", "Age", "LowerBound", "UpperBound"]
    pd.DataFrame(df).to_csv(args.nwk.replace(".nwk", ".tsv"), index=False, header=header, sep="\t")
