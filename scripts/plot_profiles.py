#!python3
import argparse
import os
import subprocess
from csv import reader
import pandas as pd
from plot_module import plot_correlation

nperline = 100
header = "site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n"
cmd_plot = "dms2_logoplot --prefs {0} --outdir {1} --name {2} --nperline {3}"
cmd_diff_plot = "dms2_logoplot --diffprefs {0} --outdir {1} --name {2} --nperline {3}"


def plot_profiles(args_input, args_infer, args_output):
    axis_dict = dict()
    if args_input != "":
        name = os.path.basename(args_input).replace(".prefs", "").replace("_", "-")
        subprocess.call(cmd_plot.format(args_input, args_output, name, nperline), shell=True)
        input_df = pd.read_csv(args_input, sep=",")
        axis_dict["Simulation"] = input_df.drop('site', axis=1).values.flatten()

    for profile in args_infer:
        prefs = args_output + "/" + os.path.basename(profile).replace(".siteprofiles", ".prefs")
        name = os.path.basename(prefs).replace(".prefs", "").replace("_", "-")
        r = reader(open(profile, 'r'), delimiter='\t')
        next(r)
        with open(prefs, 'w') as w:
            w.write(header)
            for line in r:
                w.write(",".join(line) + "\n")

        subprocess.call(cmd_plot.format(prefs, args_output, name, nperline), shell=True)
        df = pd.read_csv(prefs, sep=",")
        axis_dict[name] = df.drop('site', axis=1).values.flatten()

        if args_input != "":
            diff_df = df - input_df
            diff_df["site"] = range(1, len(df) + 1)
            diff_df.to_csv(prefs + ".diff", sep=',', encoding='utf-8', index=False)
            subprocess.call(cmd_diff_plot.format(prefs + ".diff", args_output, name + ".diff", nperline), shell=True)

    plot_correlation(os.path.join(args_output, 'correlation.aa-preferences.png'), axis_dict, {}, global_min_max=True, alpha=0.0025)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, default="", type=str, dest="input")
    parser.add_argument('-p', '--infer', required=True, type=str, nargs='+', dest="infer")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    args = parser.parse_args()
    plot_profiles(args.input, args.infer, args.output)
