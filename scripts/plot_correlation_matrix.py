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


def partial_corr_from_precision_matrix(p):
    partial_corr_matrix = np.zeros(p.shape, dtype=np.float64)

    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            partial_corr_matrix[i][j] = - p[i][j] / (np.sqrt(p[i][i] * p[j][j]))
    return partial_corr_matrix


def save_heatmap(mean_matrix, matrix_list, names, output_name):
    p_values = np.zeros(mean_matrix.shape, dtype=np.float64)
    for i in range(mean_matrix.shape[0]):
        for j in range(mean_matrix.shape[1]):
            if i == j:
                p_values[i][j] = None
            else:
                signs = [np.sign(mean_matrix[i][j]) != np.sign(m[i][j]) for m in matrix_list]
                p_values[i][j] = np.sum(signs) / len(signs)

    sns.heatmap(pd.DataFrame(p_values, index=names, columns=names), cmap='RdYlGn', annot=mean_matrix, linewidths=.5,
                vmin=0, vmax=0.05)
    plt.tight_layout()
    plt.savefig(output_name + ".svg", format='svg')
    plt.clf()
    plt.close('all')


def plot_correlation_matrix(input_trace, output_plot, burn_in):
    trace = pd.read_csv(input_trace + '.trace', sep='\t')
    if len(trace.index) <= burn_in:
        print("Burn-in was set to {0} but there is only {1} points to read.".format(burn_in, len(trace.index)))
        burn_in = max(len(trace.index) - 10, 0)
        print("Setting the burn-in to {0}.".format(burn_in))

    dim = max([int(col.split("_")[-1]) for col in trace if ("Precision_" in col or "cov_" in col)]) + 1

    precision_matrix = np.zeros((dim, dim), dtype=np.float64)
    precision_matrix_list = [np.zeros((dim, dim), dtype=np.float64) for _ in range(len(trace.index[burn_in:]))]

    for i in range(dim):
        for j in range(dim):
            if i >= j:
                name = "Precision_{0}_{1}".format(i, j)
                if name not in trace:
                    name = "cov_{0}_{1}".format(i, j)

                for point, val in enumerate(trace[name][burn_in:]):
                    precision_matrix_list[point][i][j] = val
                    if i != j:
                        precision_matrix_list[point][j][i] = val

                precision_matrix[i][j] = np.mean(trace[name][burn_in:])
                if i != j:
                    precision_matrix[j][i] = precision_matrix[i][j]

    partial_corr_matrix = partial_corr_from_precision_matrix(precision_matrix)
    partial_corr_matrix_list = [partial_corr_from_precision_matrix(p) for p in precision_matrix_list]

    corr_matrix = np.linalg.inv(precision_matrix)
    corr_matrix_list = [np.linalg.inv(p) for p in precision_matrix_list]

    if "nodeomega" in input_trace:
        names = ["LogOmega", "LogMutationRate"]
    else:
        names = ["LogNe", "LogMutationRate"]
    names += ["Maximum longevity (yrs)", "Adult weight (g)", "Female maturity (days)"]
    if dim < len(names):
        names = names[:dim]

    save_heatmap(partial_corr_matrix, partial_corr_matrix_list, names, output_plot + "_partial_correlation")
    save_heatmap(corr_matrix, corr_matrix_list, names, output_plot + "_correlation")
    save_heatmap(precision_matrix, precision_matrix_list, names, output_plot + "_precision")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, dest="trace")
    parser.add_argument('-b', '--burn_in', required=False, type=int, default=0, dest="burn_in")
    args = parser.parse_args()
    plot_correlation_matrix(args.trace, args.output, args.burn_in)
