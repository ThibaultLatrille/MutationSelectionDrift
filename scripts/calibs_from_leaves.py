#!python3
import os
import argparse
import pandas as pd
from ete3 import Tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--txt', default="../DataEmpirical/Primates/primglob.calibs.txt",
                        required=False, type=str, dest="txt")
    parser.add_argument('-t', '--tree', default="../DataEmpirical/Primates/rootedtree.nhx",
                        required=False, type=str, dest="tree")
    args = parser.parse_args()
    tree = Tree(args.tree, format=1)
    csv = pd.read_csv(args.txt, sep=' ', skiprows=[0], header=None)
    df = []
    for index, (leaf_1, leaf_2, upper, lower) in csv.iterrows():
        name = tree.get_common_ancestor([leaf_1, leaf_2]).name
        df += [[name, (lower + upper) / 2, lower, upper]]

    header = ["NodeName", "Age", "LowerBound", "UpperBound"]
    pd.DataFrame(df).to_csv(args.txt.replace(".txt", ".tsv"), index=False, header=header, sep="\t")
