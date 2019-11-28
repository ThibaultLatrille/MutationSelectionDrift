#!python3
from ete3 import Tree
import argparse
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=False, type=str,
                        default='../DataEmpirical/sp73_OrthoMam/life_history_traits.tsv', dest="tsv", metavar="<tsv>",
                        help="The tsv to be re-written")
    args = parser.parse_args()
    csv = pd.read_csv(args.tsv, sep='\t')
    df = []
    for index, row in csv.iterrows():
        line = row.tolist()
        line[0] = line[0][0].upper() + line[0].split("_")[1][:3].lower()
        df += [line]

    header = list(csv)
    pd.DataFrame(df).to_csv(args.tsv + ".abbr", index=False, header=header, sep="\t")
