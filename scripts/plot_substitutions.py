#!python3
from ete3 import Tree
from plot_module import chain, NucToColor, LineCollection
from matplotlib import rc

rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{xcolor}')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

letter_3_to_1 = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
                 'Ile': 'I', 'Pro': 'P', 'Tht': 'T', 'Phe': 'F', 'Asn': 'N',
                 'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
                 'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}

letter_3_to_colors = {'Cys': 'orange', 'Asp': 'blue', 'Ser': 'brown', 'Gln': 'chocolate', 'Lys': 'crimson',
                      'Ile': 'firebrick', 'Pro': 'forestgreen', 'Tht': 'greenyellow', 'Phe': 'indigo', 'Asn': 'hotpink',
                      'Gly': 'lavender', 'His': 'lemonchiffon', 'Leu': 'magenta', 'Arg': 'mintcream', 'Trp': 'navy',
                      'Ala': 'purple', 'Val': 'red', 'Glu': 'violet', 'Tyr': 'silver', 'Met': 'orchid'}

letter_1_to_3 = {v: k for k, v in letter_3_to_1.items()}


def plot_tree(tree, substitutions, codon_seq, aa_seq, species, outputpath, font_size=12, line_type="-",
              vt_line_width=0.5,
              hz_line_width=0.2):
    node_list = tree.iter_descendants(strategy='postorder')
    node_list = chain(node_list, [tree])

    vblankline, hblankline, = [], []

    if len(tree) < 50:
        fig = plt.figure(figsize=(16, 12))
    elif len(tree) < 70:
        fig = plt.figure(figsize=(16, 12))
    else:
        fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111)

    node_pos = dict((n2, i) for i, n2 in enumerate(tree.get_leaves()[::-1]))
    max_name_size = max(len(species[n.name]) for n in tree)
    # draw tree
    for n in node_list:
        x = sum(n2.dist for n2 in n.iter_ancestors()) + n.dist
        if n.is_leaf():
            y = node_pos[n]
            common_name = species[n.name]
            node_name = " " + common_name
            if len(common_name) != max_name_size:
                node_name += " " * (max_name_size - len(common_name))

            for i, nuc in enumerate(codon_seq[n.name]):
                ax.text(x + i * 0.008 + 0.01, y, nuc, color=NucToColor(nuc), va='center', size=font_size, zorder=50)

            for i, aa in enumerate(aa_seq[n.name]):
                ax.text(x + i * 3 * 0.008, y, letter_1_to_3[aa],
                        color=letter_3_to_colors[letter_1_to_3[aa]], va='center',
                        size=font_size, zorder=60)

            ax.text(x + 0.01, y, node_name, va='center', size=font_size, zorder=100)
        else:
            y = np.mean([node_pos[n2] for n2 in n.children])
            node_pos[n] = y

            vblankline.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
            for child in n.children:
                h = node_pos[child]
                hblankline.append(((x, h), (x + child.dist, h)))

        if n.name in substitutions:
            for _, sub in substitutions[n.name].iterrows():
                sub_x = x + sub["Time"] - n.dist
                ax.scatter(sub_x, y, s=25, marker='o', zorder=1000)
                txt = "${0} \\rightarrow {1}$".format(sub["NucFrom"], sub["NucTo"])
                ax.text(sub_x, y + 0.2, txt, va='center', ha='center', size=font_size * 0.8)
                ax.text(sub_x, y - 0.3, sub["Site"] + 1, va='center', ha='center', size=font_size * 0.8)

    ax.add_collection(LineCollection(hblankline, colors='black', linestyle=line_type, linewidth=hz_line_width * 2))
    ax.add_collection(LineCollection(vblankline, colors='black', linestyle=line_type, linewidth=vt_line_width * 2))

    # scale line
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    diffy = ymax - ymin
    dist = round((xmax - xmin) / 4, 1)
    padding = 200.
    ymin -= diffy / padding
    ax.plot([xmin, xmin + dist], [ymin, ymin], color='k')
    ax.plot([xmin, xmin], [ymin - diffy / padding, ymin + diffy / padding], color='k')
    ax.plot([xmin + dist, xmin + dist], [ymin - diffy / padding, ymin + diffy / padding],
            color='k')
    ax.text((xmin + xmin + dist) / 2, ymin - diffy / padding, dist, va='top',
            ha='center', size=font_size)

    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(outputpath, format=outputpath[outputpath.rfind('.') + 1:])
    plt.show()
    plt.close("all")


def plot_simulation(input_simu, args_output, args_species):
    t = Tree(input_simu, format=1)
    substitutions = {k: v for (k, v) in
                     pd.read_csv(input_simu.replace(".nhx", ".substitutions.tsv"), sep='\t').groupby(by="NodeName")}
    sequences = pd.read_csv(input_simu.replace(".nhx", ".sequences.tsv"), sep='\t')
    codon_seq = {row["NodeName"]: row["CodonSequence"] for _, row in sequences.iterrows()}
    aa_seq = {row["NodeName"]: row["AASequence"] for _, row in sequences.iterrows()}

    species = {row[0]: row[1] for _, row in pd.read_csv(args_species, sep='\t').iterrows()}

    plot_tree(t, substitutions, codon_seq, aa_seq, species, args_output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, default='',
                        dest="tree", metavar="<tree>", help="The tree to be drawn")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-s', '--species', required=True, type=str, dest="species")
    args = parser.parse_args()
    plot_simulation(args.tree, args.output, args.species)
