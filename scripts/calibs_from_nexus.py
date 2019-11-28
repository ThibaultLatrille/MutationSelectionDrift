#!python3
import os
import argparse
import pandas as pd
from Bio import Phylo
from ete3 import Tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--nexus', default="../DataEmpirical/Cetacea/FigTree_parts_10_mcmctree_AR.tre",
                        required=False, type=str, dest="nexus")
    parser.add_argument('-t', '--tree', default="../DataEmpirical/Cetacea/rootedtree.nhx",
                        required=False, type=str, dest="tree")
    args = parser.parse_args()
    tree = Tree(args.tree, format=1)
    nexus = Phylo.read(args.nexus, 'nexus')
    assert (len(tree) == nexus.count_terminals())
    tree_leaf_names = set(tree.get_leaf_names())
    nexus_leaf_names = set([i.name for i in nexus.get_terminals()])
    assert(nexus_leaf_names == tree_leaf_names)
    df = []
    for clade in nexus.get_nonterminals(order='postorder'):
        ages = [clade.distance(l) for l in clade.get_terminals()]
        assert(max(ages) - min(ages) < 1e-4)
        age = sum(ages) / len(ages)
        name = tree.get_common_ancestor([i.name for i in clade.get_terminals()]).name
        bounds = clade.comment[clade.comment.find("{")+1:clade.comment.find("}")].split(",")
        assert(len(bounds) == 2)
        min_age, max_age = [float(i) for i in bounds]
        assert(age > min_age)
        assert(max_age > age)
        df += [[name, age, min_age, max_age]]

    header = ["NodeName", "Age", "LowerBound", "UpperBound"]
    pd.DataFrame(df).to_csv(args.nexus.replace(".tre", ".tsv"), index=False, header=header, sep="\t")
