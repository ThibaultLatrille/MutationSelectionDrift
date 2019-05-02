#!python3
import argparse
import os
import glob
from itertools import chain
from math import floor
from ete3 import Tree
from plot_module import *
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


def round_sig(x, sig=2):
    return round(x, sig - int(floor(np.log10(abs(x)))) - 1)


def to_coord(x, y, xmin, xmax, ymin, ymax, plt_xmin, plt_ymin, plt_width, plt_height):
    x = (x - xmin) / (xmax - xmin) * plt_width + plt_xmin
    y = (y - ymin) / (ymax - ymin) * plt_height + plt_ymin
    return x, y


def format_float(x):
    if x > 0:
        return "{0:.2f}".format(x)
    else:
        return "{0:.1f}".format(x)


def plot_tree(tree, feature, outputpath, font_size=9, line_type="-", vt_line_width=0.5, hz_line_width=3.0, nbr_steps=30,
              max_circle_size=30, min_circle_size=4):
    vlinec, vlines, hlinec, hlines, nodes, nodex, nodey, ali_lines = [], [], [], [], [], [], [], []

    fig = plt.figure(figsize=(8, int(len(tree) / 2 + 4)))
    ax = fig.add_subplot(111)

    min_annot = min(float(getattr(n, feature + "_min")) for n in tree.iter_leaves() if feature in n.features)
    max_annot = max(float(getattr(n, feature + "_max")) for n in tree.iter_leaves() if feature in n.features)
    color_map = ScalarMappable(norm=Normalize(vmin=min_annot, vmax=max_annot), cmap=plt.get_cmap("gnuplot"))

    node_pos = dict((n2, i) for i, n2 in enumerate(tree.get_leaves()[::-1]))
    node_list = tree.iter_descendants(strategy='postorder')
    node_list = chain(node_list, [tree])

    max_name_size = max(len(n.name) for n in tree)
    # draw tree
    for n in node_list:
        x = sum(n2.dist for n2 in n.iter_ancestors()) + n.dist
        node_annot = float(getattr(n, feature))
        min_node_annot = float(getattr(n, feature + "_min"))
        max_node_annot = float(getattr(n, feature + "_max"))
        name = n.name
        if len(name) != max_name_size:
            name += " " * (max_name_size - len(name))
        node_name = "{0} {1} [{2},{3}]".format(name, format_float(node_annot), format_float(min_node_annot),
                                               format_float(max_node_annot))
        if n.is_leaf():
            y = node_pos[n]
            ax.text(x + 0.06, y, node_name, va='center', size=font_size)
        else:
            y = np.mean([node_pos[n2] for n2 in n.children])
            node_pos[n] = y

            # draw vertical line
            vlinec.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
            vlines.append(node_annot)

            # draw horizontal lines
            for child in n.children:
                seg_annot = node_annot
                annot_step = (float(getattr(child, feature)) - node_annot) / (nbr_steps - 1)

                h = node_pos[child]
                x_start = x
                x_step = child.dist / nbr_steps

                for seg in range(nbr_steps):
                    hlines.append(seg_annot)
                    hlinec.append(((x_start, h), (x_start + x_step, h)))
                    seg_annot += annot_step
                    x_start += x_step

        nodes.append({"min": min_node_annot, "max": max_node_annot})
        nodex.append(x)
        nodey.append(y)

    hline_col = LineCollection(hlinec, colors=[color_map.to_rgba(l) for l in hlines],
                               linestyle=line_type,
                               linewidth=hz_line_width * 2)
    vline_col = LineCollection(vlinec, colors=[color_map.to_rgba(l) for l in vlines],
                               linestyle=line_type,
                               linewidth=vt_line_width * 2)
    ali_line_col = LineCollection(ali_lines, colors='k')

    ax.add_collection(hline_col)
    ax.add_collection(vline_col)
    ax.add_collection(ali_line_col)

    def circle_radius(annot):
        return ((max_circle_size - min_circle_size) * (annot - min_annot) / (
                max_annot - min_annot) + min_circle_size) ** 2 / 2

    for zorder, suffix in enumerate(["max", "min"]):
        scat = ax.scatter(nodex, nodey, s=0, marker='o')
        scat.set_sizes([circle_radius(n[suffix]) for n in nodes])
        scat.set_color([color_map.to_rgba(n[suffix]) for n in nodes])
        scat.set_zorder(10 + zorder)

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
    color_map._A = []
    cbar = fig.colorbar(color_map, ax=ax, orientation='horizontal', pad=0, shrink=0.8)
    cbar.ax.set_xlabel(feature, labelpad=0)
    plt.tight_layout()
    plt.savefig(outputpath)
    plt.show()
    plt.close("all")


def plot_trace(input_trace, output_plot):
    axis_trees, axis_filenames = dict(), dict()

    for filepath in input_trace:
        for tree_path in sorted(glob.glob("{0}.*.nhx".format(filepath))):
            feature = tree_path.replace(filepath + ".", "").replace(".nhx", "")
            filename = os.path.basename(filepath)
            tree = Tree(tree_path, format=1)

            if feature not in axis_trees:
                axis_trees[feature] = []
            if feature not in axis_filenames:
                axis_filenames[feature] = []

            axis_filenames[feature].append(filename)
            axis_trees[feature].append(tree)
            plot_tree(tree.copy(), feature, "{0}/{1}.{2}.svg".format(output_plot, filename, feature))

    for feature in axis_trees:
        axis_dict = dict()
        for filename, tree in zip(axis_filenames[feature], axis_trees[feature]):
            axis = np.array([float(getattr(node, feature)) for node in tree.traverse()])
            if len(axis) > 0:
                axis_dict[filename] = axis

        if len(axis_dict) > 1:
            plot_correlation('{0}/correlation.{1}.png'.format(output_plot, feature), axis_dict, {}, [])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, nargs='+', dest="trace")
    args = parser.parse_args()
    plot_trace(args.trace, args.output)
