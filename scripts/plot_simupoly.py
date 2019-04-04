#!python3
from plot_module import *
import argparse
import os


def plot_simulation(input_simu):
    subs_file = input_simu.replace(".nhx", ".substitutions.tsv")
    if os.path.isfile(subs_file):
        subs = pd.read_csv(subs_file, sep='\t')
        nbr_bins = max((int(len(subs) / 100), 10))
        subs.hist(column=["selection_coefficient", "mean_fitness", "time_between", 'time_event'], bins=nbr_bins,
                  figsize=(19.20, 10.80))
        plt.savefig("{0}.png".format(subs_file))
        plt.clf()
        plt.close("all")
        subs.plot.scatter(x='time_event', y='time_between', c='selection_coefficient', figsize=(19.20, 10.80))
        plt.savefig("{0}.time.png".format(subs_file))
        plt.clf()
        plt.close("all")

    if os.path.isfile(input_simu):
        t = Tree(input_simu, format=3)
        args_nodes = set()
        for node in t.traverse():
            args_nodes = args_nodes.union(node.features)

        branch_dict = {}
        for arg in args_nodes:
            values = [to_float(getattr(n, arg)) for n in t.traverse() if
                      arg in n.features and convertible_to_float(getattr(n, arg))]
            if len(values) > 1:
                if ("Branch" in arg) and (("dNd" in arg) or ("LogNe" in arg)):
                    branch_dict[arg] = values

        der_pop_size = [(np.log10(to_float(n.population_size)) - to_float(n.Branch_LogNe)) for n in t.traverse() if
                        not n.is_root()]
        plot_correlation("{0}.correlation.png".format(input_simu), branch_dict, {}, der_pop_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, default='',
                        dest="t", metavar="<tree>", help="The tree to be drawn")
    args = parser.parse_args()
    plot_simulation(args.t)
