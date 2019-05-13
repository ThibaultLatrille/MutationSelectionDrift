#!python3
import argparse
import os
import subprocess
from csv import reader
import pandas as pd
from plot_module import plot_correlation

nperline = 100
header = "site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n"
cmd = "dms2_logoplot --prefs {0} --outdir {1} --name {2} --nperline {3}"


def plot_profiles(args_infer, args_output, axis_dict):
    for profile in args_infer:
        prefs = args_output + "/" + os.path.basename(profile).replace(".siteprofiles", ".prefs")
        name = os.path.basename(prefs).replace(".prefs", "").replace("_", "-")
        r = reader(open(profile, 'r'), delimiter='\t')
        next(r)
        with open(prefs, 'w') as w:
            w.write(header)
            for line in r:
                w.write(",".join(line) + "\n")

        subprocess.call(cmd.format(prefs, args_output, name, nperline), shell=True)
        df = pd.read_csv(prefs, sep=",")
        df.drop('site', axis=1, inplace=True)
        axis_dict[name] = df.values.flatten()

    plot_correlation(os.path.join(args_output, 'correlation.png'), axis_dict, {}, [])


def open_input(input_prefs, args_output, axis_dict):
    name = os.path.basename(input_prefs).replace(".prefs", "").replace("_", "-")
    subprocess.call(cmd.format(input_prefs, args_output, name, nperline), shell=True)
    df = pd.read_csv(input_prefs, sep=",")
    df.drop('site', axis=1, inplace=True)
    axis_dict[name] = df.values.flatten()
    return axis_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, default="", type=str, dest="input")
    parser.add_argument('-p', '--infer', required=True, type=str, nargs='+', dest="infer")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    args = parser.parse_args()
    profiles_dict = dict()
    if args.input != "":
        profiles_dict = open_input(args.input, args.output, profiles_dict)
    plot_profiles(args.infer, args.output, profiles_dict)
