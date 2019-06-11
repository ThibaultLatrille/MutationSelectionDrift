#!python3
from itertools import chain
import numpy as np
import statsmodels.api as sm
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

my_dpi = 128
RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"


def NucToColor(nuc):
    if nuc == "A":
        return RED
    elif nuc == "C":
        return YELLOW
    elif nuc == "G":
        return GREEN
    elif nuc == "T":
        return BLUE
    else:
        return LIGHTGREEN


def convertible_to_float(f):
    try:
        float(f)
        return True
    except ValueError:
        return False


def to_float(f):
    c = float(f)
    if np.isnan(c) or np.isinf(c):
        return 0.0
    else:
        return c


def format_float(x):
    if x > 0:
        return "{0:.2f}".format(x)
    else:
        return "{0:.1f}".format(x)


def tex_float(x):
    s = "{0:.3g}".format(x)
    if "e" in s:
        mantissa, exp = s.split('e')
        return mantissa + '\\times 10^{' + exp + '}'
    else:
        return s


def min_max(axis, err):
    eps = 0.05
    if isinstance(err, float):
        min_axis, max_axis = min(axis) - err, max(axis) + err
    elif err.ndim == 1:
        min_axis, max_axis = min(axis - err), max(axis + err)
    else:
        min_axis, max_axis = min(axis - err[0]), max(axis + err[1])
    min_axis -= (max_axis - min_axis) * eps
    max_axis += (max_axis - min_axis) * eps
    if min_axis == max_axis:
        return min_axis - eps, max_axis + eps
    else:
        return min_axis, max_axis


def plot_correlation(name, axis_dict, err_dict, color_map):
    f, axs = plt.subplots(len(axis_dict), len(axis_dict),
                          figsize=(len(axis_dict) * 640 / my_dpi, len(axis_dict) * 480 / my_dpi), dpi=my_dpi)

    if len(axis_dict) == 1:
        axs = [[axs]]

    scatter_kwargs = {"zorder": 0}
    error_kwargs = {"lw": .5, "zorder": -1}

    for row, (row_filename, row_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
        for col, (col_filename, col_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
            ax = axs[row][col]
            cm = plt.cm.get_cmap('Spectral_r')
            if col == 0:
                ax.set_ylabel(row_filename)
            if row == len(axis_dict) - 1:
                ax.set_xlabel(col_filename)

            min_col, max_col = min_max(col_axis, err_dict[col_filename] if (col_filename in err_dict) else 0.0)
            ax.set_xlim((min_col, max_col))

            if row == col:
                nb_bins = 50 if len(col_axis) > 100 else int(len(col_axis) / 2)
                ax.hist(col_axis, nb_bins, density=True, color=BLUE)
            else:
                min_row, max_row = min_max(row_axis, err_dict[row_filename] if (row_filename in err_dict) else 0.0)
                ax.set_ylim((min_row, max_row))

                if len(color_map) > 0 and len(color_map) == len(col_axis) and len(color_map) == len(row_axis):
                    sc = ax.scatter(col_axis, row_axis, c=color_map, cmap=cm, **scatter_kwargs,
                                    label=r"${0}$ points".format(len(row_axis)))
                else:
                    sc = ax.scatter(col_axis, row_axis, color=BLUE, **scatter_kwargs,
                                    label=r"${0}$ points".format(len(row_axis)))

                if (col_filename in err_dict) and (row_filename in err_dict):
                    ax.errorbar(col_axis, row_axis, xerr=err_dict[col_filename], yerr=err_dict[row_filename], fmt='o',
                                marker=None, mew=0, ecolor=GREEN, **error_kwargs)
                elif col_filename in err_dict:
                    ax.errorbar(col_axis, row_axis, xerr=err_dict[col_filename], fmt='o', marker=None, mew=0,
                                ecolor=GREEN, **error_kwargs)
                elif row_filename in err_dict:
                    ax.errorbar(col_axis, row_axis, yerr=err_dict[row_filename], fmt='o', marker=None, mew=0,
                                ecolor=GREEN, **error_kwargs)

                idf = np.linspace(min_col, max_col, 30)
                if row_filename != col_filename and len(set(col_axis)) > 1 and len(set(row_axis)) > 1:
                    model = sm.OLS(row_axis, sm.add_constant(col_axis))
                    results = model.fit()
                    b, a = results.params[0:2]
                    ax.plot(idf, a * idf + b, '-', color=RED, label=r"$y={0}x {3} {1}$ ($r^2={2})$".format(
                        tex_float(float(a)), tex_float(abs(float(b))), tex_float(results.rsquared),
                        "+" if float(b) > 0 else "-"))
                    if a > 0:
                        ax.plot(idf, idf, '-', color='black', label=r"$y=x$")
                    ax.legend()

    plt.tight_layout()
    if len(color_map) > 0:
        f.subplots_adjust(right=0.9)
        cbar_ax = f.add_axes([0.92, 0.05, 0.05, 0.90])
        f.colorbar(sc, cax=cbar_ax)
    plt.savefig(name, format='png')
    plt.clf()
    plt.close('all')


def plot_tree(tree, feature, outputpath, font_size=12, line_type="-", vt_line_width=0.5, hz_line_width=0.2,
              max_circle_size=20, min_circle_size=4):
    node_list = tree.iter_descendants(strategy='postorder')
    node_list = chain(node_list, [tree])

    vlinec, vlines, vblankline, hblankline, nodes, nodex, nodey = [], [], [], [], [], [], []

    if len(tree) < 50:
        fig = plt.figure(figsize=(16, 9))
    elif len(tree) < 70:
        fig = plt.figure(figsize=(16, 12))
    else:
        fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111)

    min_annot_list = [float(getattr(n, feature + "_min")) for n in tree.iter_leaves() if feature + "_min" in n.features]
    if len(min_annot_list) > 0:
        min_annot = min(min_annot_list)
    else:
        min_annot = min(float(getattr(n, feature)) for n in tree.iter_leaves() if feature in n.features)

    max_annot_list = [float(getattr(n, feature + "_max")) for n in tree.iter_leaves() if feature + "_max" in n.features]
    if len(max_annot_list) > 0:
        max_annot = max(max_annot_list)
    else:
        max_annot = max(float(getattr(n, feature)) for n in tree.iter_leaves() if feature in n.features)

    cmap = plt.get_cmap("Spectral_r")
    color_map = ScalarMappable(norm=Normalize(vmin=min_annot, vmax=max_annot), cmap=cmap)

    node_pos = dict((n2, i) for i, n2 in enumerate(tree.get_leaves()[::-1]))

    max_name_size = max(len(n.name) for n in tree)
    # draw tree
    for n in node_list:
        x = sum(n2.dist for n2 in n.iter_ancestors()) + n.dist

        min_node_annot, max_node_annot = False, False
        if (feature + "_min" in n.features) and (feature + "_max" in n.features):
            min_node_annot = float(getattr(n, feature + "_min"))
            max_node_annot = float(getattr(n, feature + "_max"))
            nodes.append({"min": min_node_annot, "max": max_node_annot})

        if n.is_leaf():
            y = node_pos[n]

            node_name = " " + n.name
            if len(n.name) != max_name_size:
                node_name += " " * (max_name_size - len(n.name))
            if feature in n.features:
                node_name += " " + format_float(float(getattr(n, feature)))
            if min_node_annot and max_node_annot:
                node_name += " [{0},{1}]".format(format_float(min_node_annot), format_float(max_node_annot))
            ax.text(x, y, node_name, va='center', size=font_size)
        else:
            y = np.mean([node_pos[n2] for n2 in n.children])
            node_pos[n] = y

            if feature in n.features:
                node_annot = float(getattr(n, feature))
                # draw vertical line
                vlinec.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
                vlines.append(node_annot)

                # draw horizontal lines
                for child in n.children:
                    child_annot = float(getattr(child, feature))
                    h = node_pos[child]
                    xs = [[x, x], [x + child.dist, x + child.dist]]
                    ys = [[h - hz_line_width, h + hz_line_width], [h - hz_line_width, h + hz_line_width]]
                    zs = [[node_annot, node_annot], [child_annot, child_annot]]
                    ax.pcolormesh(xs, ys, zs, cmap=cmap, norm=Normalize(vmin=min_annot, vmax=max_annot),
                                  shading="gouraud")
            else:
                vblankline.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
                for child in n.children:
                    h = node_pos[child]
                    hblankline.append(((x, h), (x + child.dist, h)))

        nodex.append(x)
        nodey.append(y)

    vline_col = LineCollection(vlinec, colors=[color_map.to_rgba(l) for l in vlines],
                               linestyle=line_type,
                               linewidth=vt_line_width * 2)
    ax.add_collection(LineCollection(hblankline, colors='black', linestyle=line_type, linewidth=hz_line_width * 2))
    ax.add_collection(LineCollection(vblankline, colors='black', linestyle=line_type, linewidth=vt_line_width * 2))
    ax.add_collection(vline_col)

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
    plt.savefig(outputpath, format=outputpath[-3:])
    plt.show()
    plt.close("all")
