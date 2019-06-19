#!python3
import argparse
import pandas as pd
import numpy as np
import imgkit
import matplotlib

matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

my_dpi = 128
import seaborn as sns

sns.set()


def plot_correlation_matrix(input_trace, output_plot, burn_in):
    trace = pd.read_csv(input_trace + '.trace', sep='\t')
    dim = max([int(col.split("_")[-1]) for col in trace if ("Precision_" in col or "cov_" in col)]) + 1

    precision_matrix = np.zeros((dim, dim), dtype=np.float64)

    for i in range(dim):
        for j in range(dim):
            if i >= j:
                name = "Precision_{0}_{1}".format(i, j)
                if name not in trace:
                    name = "cov_{0}_{1}".format(i, j)
                vals = trace[name]
                val = np.mean(vals[burn_in:])
                precision_matrix[i][j] = val
                if i != j:
                    precision_matrix[j][i] = val

    corr_matrix = np.linalg.inv(precision_matrix)
    for i in range(dim):
        for j in range(dim):
            corr_matrix[i][j] = - precision_matrix[i][j] / (np.sqrt(precision_matrix[i][i] * precision_matrix[j][j]))

    if "nodeomega" in input_trace:
        names = ["LogOmega", "LogMutationRate"]
    else:
        names = ["LogNe", "LogMutationRate"]
    names += ["Maximum longevity (yrs)", "Adult weight (g)", "Female maturity (days)"]
    if dim < len(names):
        names = names[:dim]
    sns.heatmap(pd.DataFrame(precision_matrix, index=names, columns=names), cmap='RdYlGn_r', annot=True, linewidths=.5)
    plt.tight_layout()
    plt.savefig(output_plot.replace(".svg", "precision.svg"), format='svg')
    plt.clf()
    plt.close('all')
    sns.heatmap(pd.DataFrame(corr_matrix, index=names, columns=names), cmap='RdYlGn_r', annot=True, linewidths=.5)
    plt.tight_layout()
    plt.savefig(output_plot, format='svg')
    plt.clf()
    plt.close('all')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, dest="trace")
    parser.add_argument('-b', '--burn_in', required=False, type=int, default=0, dest="burn_in")
    args = parser.parse_args()
    plot_correlation_matrix(args.trace, args.output, args.burn_in)
