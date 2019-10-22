#!python3
import os
import argparse
import pandas as pd
from ete3 import Tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, dest="tree")
    args = parser.parse_args()
    tree = Tree(args.tree, format=1)
    features = sorted(set.intersection(*[leaf.features for leaf in tree]), key=lambda x: x != "name")
    df = [[getattr(leaf, feature) for feature in features] for leaf in tree]
    pd.DataFrame(df).to_csv(args.tree.replace(".nhx", ".tsv"), index=False, header=features, sep="\t")
