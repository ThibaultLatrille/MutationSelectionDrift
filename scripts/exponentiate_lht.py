#!python3
import pandas as pd
import numpy as np
import argparse


def transform(x):
    try:
        return np.exp(float(x))
    except ValueError:
        return x


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, default='',
                        dest="input", metavar="<input>", help="Input file")
    args = parser.parse_args()
    traits = pd.read_csv(args.input, sep="\t")
    assert (len(traits["TaxonName"]) == len(set(traits["TaxonName"])))
    file = open(args.input.replace("traits.tsv", "exponentiate.traits.tsv"), 'w')
    file.write("#TRAITS\n{0} {1} ".format(len(traits), len(list(traits)) - 1) +
               " ".join([i.replace("Log", "") for i in list(traits)[1:]]) + "\n" +
               traits.applymap(transform).to_csv(index=False, na_rep="-1", sep=" ", header=False))
    file.close()
