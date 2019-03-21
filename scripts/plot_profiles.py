#!python3
import argparse
import os
import dms_tools2.plot
import subprocess
from csv import reader

cmd = "dms2_logoplot --prefs {0} --outdir {1} --name {2} --nperline {3}"


def plot_profiles(input_prefs, args_infer, args_output):
    nperline = 53
    header = "site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n"

    name = os.path.basename(input_prefs).replace(".prefs", "").replace("_", "-")
    subprocess.call(cmd.format(input_prefs, args_output, name, nperline), shell=True)
    prefsfiles = [input_prefs]
    names = [name]

    for profile in args_infer:
        prefs = args_output + "/" + os.path.basename(profile).replace(".siteprofiles", ".prefs")
        name = os.path.basename(prefs).replace(".prefs", "").replace("_", "-")
        r = reader(open(profile, 'r'), delimiter='\t')
        next(r)
        with open(prefs, 'w') as w:
            w.write(header)
            for line in r:
                w.write(",".join(line) + "\n")

        prefsfiles.append(prefs)
        names.append(name)
        subprocess.call(cmd.format(input_prefs, args_output, name, nperline), shell=True)

    plotfile = os.path.join(args_output, 'correlation.pdf')
    dms_tools2.plot.plotCorrMatrix(names, prefsfiles, plotfile, datatype='prefs')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input")
    parser.add_argument('-p', '--infer', required=True, type=str, nargs='+', dest="infer")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    args = parser.parse_args()
    plot_profiles(args.input, args.infer, args.output)
