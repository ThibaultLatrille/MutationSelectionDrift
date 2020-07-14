#!python3
import os
import argparse
import pandas as pd
from ete3 import Tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--nwk', default="../DataEmpirical/PrimatesBinaryLHTShort/rootedtree.nwk",
                        required=False, type=str, dest="nwk")
    args = parser.parse_args()
    nwk = Tree(args.nwk, format=1)
    root_age = nwk.get_closest_leaf()[1]
    nwk.dist = nwk.dist / root_age
    for n in nwk.iter_descendants():
        print("{0}: {1}".format(n.name, n.dist / root_age) )
        n.dist = n.dist / root_age
    nwk.write(format=1, outfile=args.nwk + ".scaled.nwk")
