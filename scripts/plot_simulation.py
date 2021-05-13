#!python3
import os
from ete3 import Tree
import argparse
from plot_module import *
from collections import defaultdict, Counter

my_dpi = 256
codon_table = defaultdict(lambda: '-')
codon_table.update({
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'})
amino_acids_set = set(codon_table.values())
amino_acids_set.remove('X')
amino_acids = "".join(sorted(amino_acids_set))
assert len(amino_acids) == 20, "There is not 20 amino-acids in the codon table"


def shanon_diversity(frequencies):
    f = [i for i in frequencies if i != 0.0]
    return -np.sum(f * np.log(f))


def open_ali_file(ali_path):
    assert (os.path.isfile(ali_path))
    species_list, ali_list = [], []
    with open(ali_path, 'r') as ali_file:
        next(ali_file)
        for line in ali_file:
            if line != "\n":
                spe, seq = line.replace("  ", " ").replace("\t", " ").replace("\n", "").split(" ")
                ali_list.append(seq)
                species_list.append(spe)
    return species_list, ali_list


def aa_freq_from_ali(input_simu):
    sp, alignment = open_ali_file(input_simu)
    alignment = np.array([list(s) for s in alignment])
    aa_frequencies = []
    for site in range(int(alignment.shape[1] / 3)):
        counter = Counter([codon_table["".join(c)] for c in alignment[:, 3 * site:3 * site + 3]])
        sum_counter = sum(counter.values()) - counter["-"]
        aa_frequencies.append([counter[aa] / sum_counter for aa in amino_acids])

    return aa_frequencies


def plot_xy(x, y, name, title):
    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    plt.scatter(x, y, c=BLUE, label="{0} points".format(len(x)))

    idf = np.linspace(min(x), max(x), 10)
    model = sm.OLS(y, sm.add_constant(x))
    results = model.fit()
    b, a = results.params[0:2]
    plt.plot(idf, a * idf + b, '-', color=RED, label=r"$y={0}x {3} {1}$ ($r^2={2})$".format(
        tex_float(float(a)), tex_float(abs(float(b))), tex_float(results.rsquared),
        "+" if float(b) > 0 else "-"))

    plt.ylim((min(y), max(y)))
    plt.legend()
    plt.xlabel('Amino-acid predicted ' + title)
    plt.ylabel("Amino-acid observed " + title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(name + ".pdf", format='pdf')
    plt.savefig(name + ".png", format='png')
    plt.clf()
    plt.close('all')


def plot_simulation(input_simu, args_output, args_prefs, args_beta):
    if args_prefs != "":
        x = pd.read_csv(args_prefs, sep=",").drop('site', axis=1).values
        for i, row in enumerate(x):
            p = np.power(row, args_beta)
            x[i, :] = p / np.sum(p)
        y = np.array(aa_freq_from_ali(input_simu.replace(".nhx", ".ali")))
        plot_xy([shanon_diversity(row) for row in x], [shanon_diversity(row) for row in y],
                "{0}/correlation.prefs.shannon".format(args_output), "entropy")

        plot_xy(x.flatten(), y.flatten(), "{0}/correlation.prefs".format(args_output), "frequency")

    t = Tree(input_simu, format=1)

    args_nodes = set()
    for node in t.traverse():
        args_nodes = args_nodes.union(node.features)

    branch_dict = {}
    for arg in args_nodes:
        if arg != "dist" and arg != "support":
            values = np.array([float(getattr(n, arg)) for n in t.traverse() if
                               arg in n.features and convertible_to_float(getattr(n, arg))])
            if len(values) > 1 and len(values) == len(list(t.traverse())):
                root_pop_size = float(getattr(t.get_tree_root(), arg))
                for n in t.traverse():
                    n.add_feature("Log" + arg, np.log(float(getattr(n, arg)) / root_pop_size))
                plot_tree(t, "Log" + arg, "{0}/tree.{1}.pdf".format(args_output, arg))
            if len(values) > 1 and ("Branch" in arg) and (("dNd" in arg) or ("LogNe" in arg)):
                branch_dict[arg] = values

    plot_correlation("{0}/correlation.Ne.dNdS.pdf".format(args_output), branch_dict, {}, global_min_max=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, default='',
                        dest="t", metavar="<tree>", help="The tree to be drawn")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-b', '--beta', required=False, type=float, default=1.0,
                        dest="beta", metavar="<beta>", help="The beta parameter")
    parser.add_argument('-p', '--prefs', required=False, type=str, default='',
                        dest="prefs", metavar="<prefs>", help="The fitness profiles")
    args = parser.parse_args()
    plot_simulation(args.t, args.output, args.prefs, args.beta)
