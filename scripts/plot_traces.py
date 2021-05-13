#!python3
import argparse
import os
from plot_module import *

my_dpi = 256


def is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def plot_trace(input_trace, output_plot):
    traces, trees = dict(), dict()
    filenames = []

    for filepath in input_trace:
        if not os.path.isfile(filepath):
            continue
        filenames.append(os.path.basename(filepath))
        df = pd.read_csv(filepath, sep='\t')
        for x_param, vals in df.items():
            if x_param not in traces:
                traces[x_param] = dict()
            traces[x_param][filenames[-1]] = vals

        branches = [c.replace("*BranchPopSize_", "") for c in df if "*BranchPopSize_" in c]
        if len(branches) == 0: continue

        data = {"Branch": [], "a": [], "r2": []}
        for branch in branches:
            burn_in = 1000
            if len(df["*BranchMutRate_" + branch].values) <= 1200:
                burn_in = 0
            x = df["*BranchMutRate_" + branch].values[burn_in:]
            y = df["*BranchPopSize_" + branch].values[burn_in:]
            plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
            plt.scatter(x, y, c=BLUE, label="{0} MCMC points".format(len(x)))

            idf = np.linspace(min(x), max(x), 10)
            model = sm.OLS(y, sm.add_constant(x))
            results = model.fit()
            if len(results.params) >= 2:
                b, a = results.params[0:2]
                plt.plot(idf, a * idf + b, '-', color=RED, label=r"$y={0}x {3} {1}$ ($r^2={2})$".format(
                    tex_float(float(a)), tex_float(abs(float(b))), tex_float(results.rsquared),
                    "+" if float(b) > 0 else "-"))

                data["Branch"].append(branch)
                data["a"].append(a)
                data["r2"].append(results.rsquared)

            plt.legend()
            plt.xlabel('Mutation rate per unit of time')
            plt.ylabel("Effective population size")
            plt.title("Branch " + branch)
            plt.legend()
            plt.tight_layout()
            plt.savefig('{0}/correlation.{1}.pdf'.format(output_plot, branch), format='pdf')
            plt.clf()
            plt.close('all')

        pd.DataFrame(data).to_csv('{0}/correlation_{1}.tsv'.format(output_plot, filenames[-1]), index=False, sep="\t")
        pd.DataFrame(data).to_latex('{0}/correlation_{1}.tex'.format(output_plot, filenames[-1]), index=False, escape=False,
                                    float_format=lambda x: ("{0:.3g}".format(x) if is_float(x) else x))
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
