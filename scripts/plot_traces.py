#!python3
import argparse
import os
import pandas as pd
import matplotlib
import numpy as np
matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
my_dpi = 128


def plot_trace(input_trace, output_plot):
    traces, trees = dict(), dict()
    filenames = []

    for filepath in input_trace:
        if not os.path.isfile(filepath):
            continue
        filenames.append(os.path.basename(filepath))
        for x_param, vals in pd.read_csv(filepath, sep='\t').items():
            if x_param not in traces:
                traces[x_param] = dict()
            traces[x_param][os.path.basename(filepath)] = vals

    mean_dict = {"Name": filenames}
    for p, v in traces.items():
        mean_dict[p] = [np.mean(v[f]) if f in v else "NaN" for f in filenames]
        mean_dict[p + "Var"] = [1.96 * np.std(v[f]) if f in v else "NaN" for f in filenames]
    pd.DataFrame(mean_dict).to_csv('{0}/traces.tsv'.format(output_plot), index=False, header=mean_dict.keys(), sep="\t")

    for x_param, x_traces_param in traces.items():
        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        for x_filename, x_param_trace in sorted(x_traces_param.items(), key=lambda x: x[0]):
            style = "-"
            if "False" in x_filename:
                style = "--"
            plt.plot(range(len(x_param_trace)), x_param_trace, style, alpha=0.5, linewidth=1, label=x_filename)

        plt.xlabel('Point')
        plt.ylabel(x_param)
        plt.legend()
        plt.tight_layout()
        plt.savefig('{0}/trace.{1}.svg'.format(output_plot, x_param), format='svg')
        plt.clf()
        plt.close('all')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, nargs='+', dest="trace")
    args = parser.parse_args()
    plot_trace(args.trace, args.output)
