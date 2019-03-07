#!python3
import argparse
from ete3 import Tree, TreeStyle, TextFace, faces, CircleFace
import numpy as np
import pandas as pd
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


def layout(node, arg, min_arg, max_arg):
    if arg in node.features:
        if min_arg == max_arg:
            radius = 15
        else:
            radius = 15 * (getattr(node, arg) - min_arg) / (max_arg - min_arg) + 5
        circle = CircleFace(radius=radius, color="RoyalBlue", style="sphere")
        circle.opacity = 0.3
        faces.add_face_to_node(circle, node, 0, position="float")
        faces.add_face_to_node(TextFace(str(getattr(node, arg)) + " "), node, 0, position="aligned")


def tree_plot(input_simu):
    t = Tree(input_simu, format=3)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.complete_branch_lines_when_necessary = False

    nodes = pd.read_csv(input_simu.replace(".nhx", ".nodes.tsv"), sep='\t')
    ss_nodes = ["dN", "dS", "dNdS"]
    leaves = pd.read_csv(input_simu.replace(".nhx", ".leaves.tsv"), sep='\t')
    ss_leaves = []
    for method in ["pairwise", "watterson"]:
        for subset in ["", "_syn", "_non_syn"]:
            ss_leaves.append("theta_{0}{1}".format(method, subset))
    data = {"nodes": (nodes, ss_nodes), "leaves": (leaves, ss_leaves)}

    for csv, ss_list in data.values():
        for ss in ss_list:
            cols, groups = zip(*csv.groupby(by="taxon_name"))
            values = np.array([g[ss].values for g in groups]).T
            df = pd.DataFrame(values, columns=cols)
            df.boxplot(column=list(cols))
            for key in [ss, ss + "_pred"]:
                values = [float(getattr(n, key)) for n in t.traverse() if key in n.features]
                if len(values) > 0 and len(values) == len(cols):
                    plt.scatter(cols, values, label="Exome " + key)

            if "theta" in ss:
                plt.yscale("log")
            plt.ylabel(ss)
            plt.xticks(rotation='vertical')
            plt.legend()
            plt.savefig("{0}.tsv.{1}.png".format(input_simu.replace(".nhx", ""), ss), bbox_inches='tight')
            plt.clf()
            plt.close("all")

    args_nodes = set()
    for n in t.traverse():
        args_nodes = args_nodes.union(n.features)
    print(args_nodes)

    for arg in args_nodes:
        values = [getattr(n, arg) for n in t.traverse() if
                  arg in n.features and isinstance(getattr(n, arg), (int, float))]
        if len(values) > 0:
            max_arg = max(values)
            min_arg = min(values)
            ts.layout_fn = lambda x: layout(x, arg, min_arg, max_arg)
            for col, name in enumerate(args_nodes):
                nameF = TextFace(name, fsize=7)
                nameF.rotation = -90
                ts.aligned_header.add_face(nameF, column=col)
            ts.title.add_face(TextFace("{0} in simulation".format(arg), fsize=20), column=0)
            t.render("{0}.{1}.png".format(input_simu, arg), tree_style=ts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, default='',
                        dest="t", metavar="<tree>", help="The tree to be drawn")
    args = parser.parse_args()
    tree_plot(args.t)
