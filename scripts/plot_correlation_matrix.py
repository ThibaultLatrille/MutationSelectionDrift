#!python3
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib

matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from subprocess import run
my_dpi = 128
# 'RdYlGn'
color_map = 'RdYlGn'
sns.set()


def partial_cov_from_precision_matrix(p):
    partial_cov_matrix = np.zeros(p.shape, dtype=np.float64)
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            partial_cov_matrix[i][j] = - p[i][j] / (np.sqrt(p[i][i] * p[j][j]))
    return partial_cov_matrix


def corr_from_cov_matrix(p):
    corr_matrix = np.zeros(p.shape, dtype=np.float64)
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            corr_matrix[i][j] = p[i][j] / (np.sqrt(p[i][i] * p[j][j]))
    return corr_matrix


def save_heatmap(mean_matrix, matrix_list, names, output_name):
    p_values = np.zeros(mean_matrix.shape, dtype=np.float64)
    for i in range(mean_matrix.shape[0]):
        for j in range(mean_matrix.shape[1]):
            if i == j:
                p_values[i][j] = None
            else:
                signs = np.sum([np.sign(mean_matrix[i][j]) != np.sign(m[i][j]) for m in matrix_list])
                p_values[i][j] = signs / len(matrix_list)

    sns.heatmap(pd.DataFrame(p_values, index=names, columns=names), cmap=color_map, annot=mean_matrix, linewidths=.5,
                vmin=0.001, vmax=0.05)
    plt.savefig(output_name + ".pvalues.svg", format='svg')
    plt.clf()
    if np.max(np.abs(mean_matrix)) <= 1:
        for i in range(mean_matrix.shape[0]):
            mean_matrix[i][i] = None
        sns.heatmap(pd.DataFrame(mean_matrix, index=names, columns=names), annot=mean_matrix, cmap=color_map,
                    linewidths=.5, vmin=-1, vmax=1)
    else:
        sns.heatmap(pd.DataFrame(mean_matrix, index=names, columns=names), annot=mean_matrix, cmap=color_map,
                    linewidths=.5)
    plt.savefig(output_name + ".svg", format='svg')
    plt.clf()
    plt.close('all')


def tex_f(f):
    if 1e-3 < abs(f) < 1e3:
        return "{0:.2g}".format(f)
    else:
        base, exponent = "{0:.2e}".format(f).split("e")
        return r"{0} \times 10^{{{1}}}".format("{0:.2g}".format(float(base)), int(exponent))


def format_header(s):
    if s == "LogNe":
        return "$N_{\\mathrm{e}}$"
    elif s == "LogOmega":
        return "$\\omega$"
    elif s == "LogMutationRate":
        return "$\\mu$"
    else:
        return "\\\\".join(s.split("_"))


def save_latex(mean_matrix, matrix_list, names, output_name):
    table = open(output_name + ".tex", 'w')
    table.writelines("\\documentclass[USLetter,5pt]{article}\n"
                     "\\usepackage{adjustbox}\n")
    table.writelines("\\newcommand{\\specialcell}[2][c]{%\n\\begin{tabular}[#1]{@{}c@{}}#2\\end{tabular}}\n")
    table.writelines("\\begin{document}\n")

    heading = [""] + ["\\specialcell{" + format_header(n) + "}" for n in names]
    table.writelines("\\begin{table}[ht]\n\\centering\n\\begin{adjustbox}{width = 1\\textwidth}\n")
    table.writelines("\\begin{tabular}{|c|" + "c" * len(names) + "|}\n")
    table.writelines("\\hline\n")
    table.writelines(" & ".join(heading) + "\\\\\n")
    table.writelines("\\hline\n")
    threshold_1 = 0.05
    threshold_2 = 0.025
    assert (threshold_2 < threshold_1)
    for i in range(mean_matrix.shape[0]):
        elts = ["\\specialcell{" + format_header(names[i]) + "}"]
        for j in range(mean_matrix.shape[1]):
            if i != j:
                pval = np.sum([np.sign(mean_matrix[i][j]) != np.sign(m[i][j]) for m in matrix_list]) / len(matrix_list)
                elt = "$" + tex_f(mean_matrix[i][j])
                if pval < threshold_1:
                    elt += "^{*"
                    if pval < threshold_2:
                        elt += "*"
                    elt += "}"
                elts.append(elt + "$")
            else:
                elts.append("$\\dots$")
        table.writelines(" & ".join(elts) + "\\\\\n")
    table.writelines("\\hline\n")
    table.writelines("\\end{tabular}\n")
    table.writelines("\\end{adjustbox}\n" +
                     "\\caption{" +
                     "Asterisks indicate strength of support ($^{*} pp > " + "{0}".format(1 - threshold_1) +
                     "$, $^{**} pp > " + "{0}".format(1 - threshold_2) +
                     "$)}\n" + "\\end{table}\n")
    table.writelines("\\end{document}\n")
    table.close()
    run("pdflatex -output-directory {0} {1}.tex".format(os.path.dirname(output_name), output_name), shell=True)


def plot_covariance_matrix(input_trace, output_plot, burn_in):
    trace = pd.read_csv(input_trace + '.trace', sep='\t')
    if len(trace.index) <= burn_in:
        print("Burn-in was set to {0} but there is only {1} points to read.".format(burn_in, len(trace.index)))
        burn_in = max(len(trace.index) - 10, 0)
        print("Setting the burn-in to {0}.".format(burn_in))

    dim = max([int(col.split("_")[-1]) for col in trace if ("Precision_" in col or "cov_" in col)]) + 1

    lht_header = []
    if os.path.isfile('{0}/life_history_traits.tsv'.format(os.path.dirname(input_trace))):
        df = pd.read_csv('{0}/life_history_traits.tsv'.format(os.path.dirname(input_trace)), sep='\t')
        lht_header = list(df.columns.values)[1:]

    assert (len(lht_header) + 2 == dim)
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

    cov_matrix = np.linalg.inv(precision_matrix)
    cov_matrix_list = [np.linalg.inv(p) for p in precision_matrix_list]

    corr_matrix = corr_from_cov_matrix(cov_matrix)
    corr_matrix_list = [corr_from_cov_matrix(c) for c in cov_matrix_list]

    if "nodeomega" in input_trace:
        names = ["LogOmega", "LogMutationRate"]
    else:
        names = ["LogNe", "LogMutationRate"]
    names += lht_header

    save_latex(corr_matrix, corr_matrix_list, names, output_plot + "_correlation")

    save_heatmap(cov_matrix, cov_matrix_list, names, output_plot + "_covariance")
    save_heatmap(corr_matrix, corr_matrix_list, names, output_plot + "_correlation")
    save_heatmap(precision_matrix, precision_matrix_list, names, output_plot + "_precision")
    for in_dim in range(2, dim + 1):
        sub_par_matrix = partial_cov_from_precision_matrix(np.linalg.inv(cov_matrix[:in_dim, :in_dim]))
        sub_par_matrix_list = [partial_cov_from_precision_matrix(np.linalg.inv(c[:in_dim, :in_dim])) for c in
                               cov_matrix_list]
        save_heatmap(sub_par_matrix, sub_par_matrix_list, names[:in_dim],
                     output_plot + "_partial_correlation_{0}".format(in_dim))
        save_latex(sub_par_matrix, sub_par_matrix_list, names[:in_dim],
                   output_plot + "_partial_correlation_{0}".format(in_dim))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, dest="trace")
    parser.add_argument('-b', '--burn_in', required=False, type=int, default=0, dest="burn_in")
    args = parser.parse_args()
    plot_covariance_matrix(args.trace, args.output, args.burn_in)
