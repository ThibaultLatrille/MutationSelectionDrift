#!python3
from ete3 import Tree, TreeStyle, TextFace, faces, CircleFace
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

my_dpi = 128
RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"


def tex_f(f):
    if 0.1 < f < 100:
        return "{0:.2g}".format(f)
    else:
        return "{0:.2e}".format(f)


def layout_circle(node, arg, min_arg, max_arg):
    if min_arg == max_arg:
        radius = 15
    else:
        radius = 15 * (getattr(node, arg) - min_arg) / (max_arg - min_arg) + 5
    circle = CircleFace(radius=radius, color="RoyalBlue", style="sphere")
    circle.opacity = 0.3
    faces.add_face_to_node(circle, node, 0, position="float")


def layout(node, arg, min_arg, max_arg):
    if arg in node.features:
        layout_circle(node, arg, min_arg, max_arg)
        faces.add_face_to_node(TextFace(tex_f(getattr(node, arg)) + " "), node, 0, position="aligned")


def mutiple_layout(node, arg, min_arg, max_arg, filename, columns):
    feature = arg + "." + filename
    if feature in node.features:
        layout_circle(node, feature, min_arg, max_arg)
        for col, attr in enumerate(columns):
            faces.add_face_to_node(TextFace(tex_f(float(getattr(node, arg + "." + attr))) + " "), node, col,
                                   position="aligned")


def is_ultrametric(tree):
    d = [tree.get_distance(leaf) for leaf in tree.iter_leaves()]
    if abs(max(d) - min(d)) / max(d) < 1e-4:
        print("The tree is ultrametric.")
        return True
    else:
        print("The tree is not ultrametric.")
        print("Max distance from root to leaves: {0}.".format(max(d)))
        print("Min distance from root to leaves: {0}.".format(min(d)))
        return False


def min_max(axis, err):
    eps = 0.05
    min_axis, max_axis = min(axis - err), max(axis + err)
    min_axis -= (max_axis - min_axis) * eps
    max_axis += (max_axis - min_axis) * eps
    if min_axis == max_axis:
        return min_axis - eps, max_axis + eps
    else:
        return min_axis, max_axis


def plot_correlation(name, axis_dict, err_dict):
    f, axs = plt.subplots(len(axis_dict), len(axis_dict),
                          figsize=(len(axis_dict) * 640 / my_dpi, len(axis_dict) * 480 / my_dpi), dpi=my_dpi)

    if len(axis_dict) == 1:
        axs = [[axs]]

    for key, val in axis_dict.items():
        if key not in err_dict:
            err_dict[key] = np.zeros(len(val))

    for row, (row_filename, row_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
        for col, (col_filename, col_axis) in enumerate(sorted(axis_dict.items(), key=lambda x: x[0])):
            ax = axs[row][col]

            ax.errorbar(col_axis, row_axis,
                        xerr=err_dict[col_filename], yerr=err_dict[row_filename],
                        fmt='o', color=BLUE, ecolor=GREEN, label=r"${0}$ nodes".format(len(row_axis)))

            min_col, max_col = min_max(col_axis, err_dict[col_filename])
            idf = np.linspace(min_col, max_col, 30)
            if col == 0:
                ax.set_ylabel(row_filename)
            if row == len(axis_dict) - 1:
                ax.set_xlabel(col_filename)
            ax.set_xlim((min_col, max_col))
            min_row, max_row = min_max(row_axis, err_dict[row_filename])
            ax.set_ylim((min_row, max_row))
            if row_filename != col_filename and len(set(col_axis)) > 1 and len(set(row_axis)) > 1:
                model = sm.OLS(row_axis, sm.add_constant(col_axis))
                results = model.fit()
                b, a = results.params[0:2]
                ax.plot(idf, a * idf + b, '-', color=RED, label=r"$y={0:.3g}x {3} {1:.3g}$ ($r^2={2:.3g})$".format(
                    float(a), abs(float(b)), results.rsquared, "+" if float(b) > 0 else "-"))
                ax.legend()
    plt.tight_layout()
    plt.savefig(name, format='png')
    plt.clf()
    plt.close('all')
