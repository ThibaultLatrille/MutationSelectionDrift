#!python3
from itertools import chain
import numpy as np
import statsmodels.api as sm
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable
import matplotlib.colors as colors

matplotlib.rcParams["font.family"] = ["Latin Modern Sans"]
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.font_manager as font_manager

mono_font = font_manager.FontProperties(family='Latin Modern Mono', style='normal')
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
    if 0.001 < x < 10:
        return "{:6.3f}".format(x)
    elif x < 10000:
        return "{:6.1f}".format(x)
    else:
        s = "{:6.2g}".format(x)
        if "e" in s:
            mantissa, exp = s.split('e')
            s = mantissa + 'e$^{' + str(int(exp)) + '}$'
            s = " " * (5 + 6 - len(s)) + s
        return s


def tex_float(x):
    s = "{0:.3g}".format(x)
    if "e" in s:
        mantissa, exp = s.split('e')
        return mantissa + 'e^{' + exp + '}'
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
    # min_axis -= (max_axis - min_axis) * eps
    # max_axis += (max_axis - min_axis) * eps
    if min_axis == max_axis:
        return min_axis - eps, max_axis + eps
    else:
        return min_axis, max_axis


def label_transform(s):
    if s == "LogPopulationSize" or s == "LogNe":
        return 'Effective population size ($N_{\\mathrm{e}} $)'
    elif s == "LogMutationRate" or s == "LogMutationRatePerTime":
        return 'Mutation rate per unit of time ($\\mu$)'
    elif s == "TraitsLogGenomeSize":
        return "Genome size (Mb)"
    elif s == "LogOmega":
        return 'Non-synonymous relative substitution rate ($\\omega$)'
    elif s == "TraitsAdult_weight_":
        return 'Adult weight (g)'
    elif s == "TraitsFemale_maturity_":
        return 'Female maturity (days)'
    elif s == "TraitsMaximum_longevity_":
        return 'Maximum longevity (years)'
    elif s == "Traitsgeneration_time":
        return 'Generation time (days)'
    elif s == "Traitslongevity":
        return 'Longevity (years)'
    elif s == "Traitsmass":
        return 'Mass (kg)'
    elif s == "Traitsmaturity":
        return 'Maturity (days)'
    elif s == "TraitspiNpiS":
        return '$\\pi_{N} / \\pi_{S}$'
    elif s == "TraitspiS":
        return '$\\pi_{S}$'
    elif s == "aa-preferences":
        return 'Amino-acid preferences'
    elif s == "ContrastPopulationSize":
        return 'Contrast Population Size'
    elif s == "Log10BranchLength":
        return 'Branch length'
    elif s == "BranchTime":
        return 'Branch time'
    else:
        return s


def label_id(s):
    if "Id" in s:
        return "{0}".format(int(s.replace("Id", "")) + 1)
    if "run" in s:
        return s.split("_")[-2]
    else:
        return s


def label_corr_transform(s):
    if "read" in s or "run" in s:
        s = "Chain " + label_id(s)
    elif "Id" in s:
        s = "Replicate " + label_id(s)
    return s


def expo(x, base10=False):
    if base10:
        return np.power(10, x)
    else:
        return np.exp(x)


def plot_correlation(name, axis_dict, err_dict, global_min_max=False, alpha=1.0, font_size=14):
    if len(axis_dict) == 1: return

    is_log = "Log" in name.split(".")[-2] or "Traits" in name.split(".")[-2]
    base10 = "Log10" in name.split(".")[-2]
    scatter_kwargs = {"zorder": 0}
    error_kwargs = {"lw": .5, "zorder": -1}

    if is_log:
        for axis in axis_dict:
            axis_dict[axis] = expo(axis_dict[axis], base10)
            if axis in err_dict:
                if err_dict[axis].ndim == 1:
                    err_dict[axis] = axis_dict[axis] * (expo(err_dict[axis], base10) - 1)
                else:
                    err_dict[axis][0] = axis_dict[axis] * (1 - expo(-err_dict[axis][0], base10))
                    err_dict[axis][1] = axis_dict[axis] * (expo(err_dict[axis][1], base10) - 1)

    min_max_axis = {n: min_max(a, err_dict[n] if (n in err_dict) else 0.0) for n, a in axis_dict.items()}
    if global_min_max:
        global_min = min([min(i) for i in min_max_axis.values()])
        global_max = max([max(i) for i in min_max_axis.values()])
        for n in min_max_axis:
            min_max_axis[n] = (global_min, global_max)

    for row, (row_filename, row_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
        for col, (col_filename, col_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if row_filename == "Simulation": continue
            if col_filename != "Simulation" and col_filename >= row_filename: continue
            ax.set_xlim(min_max_axis[col_filename])
            ax.set_ylim(min_max_axis[row_filename])

            ax.scatter(col_axis, row_axis, color=BLUE, **scatter_kwargs,
                       label=r"${0}$ points".format(len(row_axis)), alpha=alpha)

            if (col_filename in err_dict) and (row_filename in err_dict):
                ax.errorbar(col_axis, row_axis, xerr=err_dict[col_filename], yerr=err_dict[row_filename], fmt='o',
                            marker=None, mew=0, ecolor=GREEN, **error_kwargs)
            elif col_filename in err_dict:
                ax.errorbar(col_axis, row_axis, xerr=err_dict[col_filename], fmt='o', marker=None, mew=0,
                            ecolor=GREEN, **error_kwargs)
            elif row_filename in err_dict:
                ax.errorbar(col_axis, row_axis, yerr=err_dict[row_filename], fmt='o', marker=None, mew=0,
                            ecolor=GREEN, **error_kwargs)
            if is_log:
                idf = np.logspace(np.log10(min(min_max_axis[col_filename])),
                                  np.log10(max(min_max_axis[col_filename])), 30)
                if row_filename != col_filename and len(set(col_axis)) > 1 and len(set(row_axis)) > 1:
                    model = sm.OLS(np.log(row_axis), sm.add_constant(np.log(col_axis)))
                    results = model.fit()
                    b, a = results.params[0:2]
                    ax.plot(idf, np.exp(a * np.log(idf) + b), '-', color=RED,
                            label=r"$y={0}x {3} {1}$ ($r^2={2})$".format(
                                tex_float(float(a)), tex_float(abs(float(b))), tex_float(results.rsquared),
                                "+" if float(b) > 0 else "-"))
                    if a > 0:
                        ax.plot(idf, idf, '-', color='black', label=r"$y=x$")
                    ax.legend(fontsize=font_size)
            else:
                idf = np.linspace(min(min_max_axis[col_filename]), max(min_max_axis[col_filename]), 30)
                if row_filename != col_filename and len(set(col_axis)) > 1 and len(set(row_axis)) > 1:
                    model = sm.OLS(row_axis, sm.add_constant(col_axis))
                    results = model.fit()
                    b, a = results.params[0:2]
                    ax.plot(idf, a * idf + b, '-', color=RED, label=r"$y={0}x {3} {1}$ ($r^2={2})$".format(
                        tex_float(float(a)), tex_float(abs(float(b))), tex_float(results.rsquared),
                        "+" if float(b) > 0 else "-"))
                    if a > 0:
                        ax.plot(idf, idf, '-', color='black', label=r"$y=x$")
                    ax.legend(fontsize=font_size)

            label = label_transform(name.split(".")[-2])
            ax.set_ylabel(label + " - " + label_corr_transform(row_filename), size=font_size * 1.2)
            ax.set_xlabel(label + " - " + label_corr_transform(col_filename), size=font_size * 1.2)
            ax.tick_params(axis='both', which='major', labelsize=font_size * 1.2)
            ax.tick_params(axis='both', which='minor', labelsize=font_size)
            if is_log:
                ax.set_xscale("log")
                ax.set_yscale("log")
            dot = name.rfind('.')
            plt.tight_layout()
            plt.savefig("{}-{}-{}.{}".format(name[:dot],
                                             col_filename if col_filename == "Simulation" else label_id(col_filename),
                                             label_id(row_filename), name[dot + 1:]), format=name[dot + 1:])
            plt.clf()
            plt.close('all')


def get_annot(n, f):
    return np.exp(float(getattr(n, f)))


def plot_tree(tree, feature, outputpath, font_size=14, line_type="-", vt_line_width=0.5, hz_line_width=0.2,
              max_circle_size=20, min_circle_size=4):
    node_list = tree.iter_descendants(strategy='postorder')
    node_list = chain(node_list, [tree])
    print(feature + "\n" + outputpath)
    vlinec, vlines, vblankline, hblankline, nodes, nodex, nodey = [], [], [], [], [], [], []

    if len(tree) < 50:
        fig = plt.figure(figsize=(16, 9))
    elif len(tree) < 70:
        fig = plt.figure(figsize=(16, 12))
    else:
        fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111)

    min_annot = min(get_annot(n, feature) for n in tree.iter_leaves() if feature in n.features)
    max_annot = max(get_annot(n, feature) for n in tree.iter_leaves() if feature in n.features)

    cmap = plt.get_cmap("inferno")
    norm = colors.LogNorm(vmin=min_annot, vmax=max_annot)
    color_map = ScalarMappable(norm=norm, cmap=cmap)

    node_pos = dict((n2, i) for i, n2 in enumerate(tree.get_leaves()[::-1]))

    max_name_size = max(len(n.name) for n in tree)
    # draw tree
    rows = []
    for n in node_list:
        x = sum(n2.dist for n2 in n.iter_ancestors()) + n.dist

        min_node_annot, max_node_annot = False, False
        if (feature + "_min" in n.features) and (feature + "_max" in n.features):
            min_node_annot = get_annot(n, feature + "_min")
            max_node_annot = get_annot(n, feature + "_max")

        if n.is_leaf():
            y = node_pos[n]

            node_name = " " + n.name
            row = {"Taxon": n.name}
            if len(n.name) != max_name_size:
                node_name += " " * (max_name_size - len(n.name))
            if feature in n.features:
                node_name += " " + format_float(get_annot(n, feature))
                row[feature] = get_annot(n, feature)
            if min_node_annot and max_node_annot:
                node_name += " [{0},{1}]".format(format_float(min_node_annot), format_float(max_node_annot))
                row[feature + "Lower"] = min_node_annot
                row[feature + "Upper"] = max_node_annot
            ax.text(x, y, node_name, va='center', size=font_size, name="Latin Modern Mono")
            rows.append(row)
        else:
            y = np.mean([node_pos[n2] for n2 in n.children])
            node_pos[n] = y

            if feature in n.features:
                node_annot = get_annot(n, feature)
                # draw vertical line
                vlinec.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
                vlines.append(node_annot)

                # draw horizontal lines
                for child in n.children:
                    child_annot = get_annot(child, feature)
                    h = node_pos[child]
                    xs = [[x, x], [x + child.dist, x + child.dist]]
                    ys = [[h - hz_line_width, h + hz_line_width], [h - hz_line_width, h + hz_line_width]]
                    zs = [[node_annot, node_annot], [child_annot, child_annot]]
                    ax.pcolormesh(xs, ys, zs, cmap=cmap, norm=norm, shading="gouraud")
            else:
                vblankline.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
                for child in n.children:
                    h = node_pos[child]
                    hblankline.append(((x, h), (x + child.dist, h)))

        nodex.append(x)
        nodey.append(y)

    pd.DataFrame(rows).to_csv(outputpath[:-4] + ".tsv", index=None, header=rows[0].keys(), sep='\t')
    vline_col = LineCollection(vlinec, colors=[color_map.to_rgba(l) for l in vlines],
                               linestyle=line_type,
                               linewidth=vt_line_width * 2)
    ax.add_collection(LineCollection(hblankline, colors='black', linestyle=line_type, linewidth=hz_line_width * 2))
    ax.add_collection(LineCollection(vblankline, colors='black', linestyle=line_type, linewidth=vt_line_width * 2))
    ax.add_collection(vline_col)

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
    # plt.tight_layout()
    color_map._A = []
    ticks = "PopulationSize" in feature or "MutationRatePerTime" in feature or "Omega" in feature
    cbar = fig.colorbar(color_map, ticks=(
        [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50] if ticks else None),
                        orientation='horizontal', pad=0, shrink=0.6)
    cbar.ax.xaxis.set_tick_params('major', labelsize=font_size * 1.8)
    cbar.ax.xaxis.set_tick_params('minor', labelsize=0)
    if ticks:
        cbar.ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    cbar.ax.set_xlabel(label_transform(feature), labelpad=5, size=font_size * 1.8)
    plt.tight_layout()
    for o in fig.findobj():
        o.set_clip_on(False)
    plt.savefig(outputpath, format=outputpath[outputpath.rfind('.') + 1:], bbox_inches='tight')
    plt.close("all")
