#!python3
import argparse
import pandas as pd
import numpy as np
import os
from subprocess import run


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


def tex_f(f):
    return "{0:.3g}".format(f)


def format_header(s):
    if s == "LogNe":
        return "$\\bm{N_{\\mathrm{e}}}$"
    elif s == "LogOmega":
        return "$\\bm{\\omega}$"
    elif s == "LogMutationRate":
        return "$\\bm{\\mu}$"
    elif s == "piS":
        return "$\\bm{\\pi_{S}}$"
    elif s == "piNpiS":
        return "$\\bm{\\pi_{N}/\\pi_{S}}$"
    else:
        if ('(' in s) and ("$" not in s): s = s.split('(')[0]
        return "\\textbf{" + s.replace('_', " ") + "}"


def save_latex(output_prefix, matrix_list, names, output_name, plot_diag=False, disp_pval=True):
    tabu = "\\centering\n\\noindent\\adjustbox{max width=\\textwidth}{%\n"
    heading = [format_header(n) for n in [output_name] + names]
    tabu += "\\begin{tabu}{|c||" + "c|" * len(names) + "}\n"
    tabu += "\\hline\n"
    tabu += " & ".join(heading) + "\\\\\n"
    tabu += "\\hhline{|=#" + "=|" * len(names) + "}\n"
    threshold_1 = 0.05
    threshold_2 = 0.025
    assert (threshold_2 < threshold_1)

    for i in range(matrix_list[0].shape[0]):
        elts = [format_header(names[i])]
        for j in range(matrix_list[0].shape[1]):
            if i < j or (i == j and plot_diag):
                x = np.mean([m[i][j] for m in matrix_list])
                pval = np.sum([np.sign(x) != np.sign(m[i][j]) for m in matrix_list]) / len(matrix_list)
                elt = "$" + tex_f(x)
                if disp_pval and pval < threshold_1:
                    elt += "^{*"
                    if pval < threshold_2:
                        elt += "*"
                    elt += "}"
                elts.append(elt + "$")
            else:
                elts.append("-")
        tabu += " & ".join(elts) + "\\\\\\hline\n"
    tabu += "\\end{tabu}}\n"
    if '(' in output_name: output_name = output_name.split('(')[0]
    f = open(output_prefix + "_" + output_name.strip().lower().replace(" ", "_") + ".tex", 'w')
    f.writelines(tabu)
    f.close()
    table = "\\begin{table}[ht]\n" + tabu
    table += "\\caption{" + output_name + ". "
    if disp_pval:
        table += "Asterisks indicate strength of support ($^{*} pp > " + "{0}".format(
            1 - threshold_1) + "$, $^{**} pp > " + "{0}$)".format(1 - threshold_2)
    table += "}\n\\end{table}\n"
    return table


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

    if len(lht_header) + 2 == dim:
        if "nodeomega" in input_trace:
            names = ["LogOmega", "LogMutationRate"]
        else:
            names = ["LogNe", "LogMutationRate"]
        names += lht_header
    elif dim == 3:
        names = ["LogNe", "LogMutationRate", "LogGenerationTime"]
    elif dim == 6:
        names = ["maturity", "mass", "longevity", "piS", "piNpiS", "generation_time"]
    elif dim == 4:
        names = ["maturity", "mass", "longevity", "generation_time"]
    else:
        exit(1)

    table = open(output_plot + ".tex", 'w')
    table.writelines("\\documentclass[USLetter,5pt]{article}\n"
                     "\\usepackage{tabu}\n"
                     "\\usepackage{bm}\n"
                     "\\usepackage{hhline}\n"
                     "\\usepackage{adjustbox}\n")
    table.writelines("\\begin{document}\n")

    precision_matrix_list = [np.zeros((dim, dim), dtype=np.float64) for _ in range(len(trace.index[burn_in:]))]
    for i in range(dim):
        for j in range(dim):
            if i >= j:
                name = "Precision_{0}_{1}".format(i, j)
                if name not in trace:
                    name = "cov_{0}_{1}".format(i, j)

                for point, val in enumerate(trace[name][burn_in:]):
                    precision_matrix_list[point][i][j] = val
                    precision_matrix_list[point][j][i] = val

    print(names)
    print(dim)
    cov_matrix_list = [np.linalg.inv(p) for p in precision_matrix_list]
    table.write(save_latex(output_plot, cov_matrix_list, names, "Covariance ($\\bm{\\Sigma}$)", True))

    corr_matrix_list = [corr_from_cov_matrix(c) for c in cov_matrix_list]
    table.write(save_latex(output_plot, corr_matrix_list, names, "Correlation ($\\bm{\\rho}$)"))

    r2_matrix_list = [np.multiply(c, c) for c in corr_matrix_list]
    table.write(save_latex(output_plot, r2_matrix_list, names, "R-squared ($\\bm{r^2}$)", False, False))

    table.write(save_latex(output_plot, precision_matrix_list, names, "Precision ($\\bm{\\Omega}$)", True))

    sub_par_matrix_list = [partial_cov_from_precision_matrix(np.linalg.inv(c)) for c in cov_matrix_list]
    table.write(save_latex(output_plot, sub_par_matrix_list, names, "Partial coefficient"))

    table.writelines("\\end{document}\n")
    table.close()
    run("pdflatex -output-directory {0} {1}.tex".format(os.path.dirname(output_plot), output_plot), shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, dest="trace")
    parser.add_argument('-b', '--burn_in', required=False, type=int, default=0, dest="burn_in")
    args = parser.parse_args()
    plot_covariance_matrix(args.trace, args.output, args.burn_in)
