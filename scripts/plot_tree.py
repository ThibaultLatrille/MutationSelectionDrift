#!python3

from ete3 import Tree, TreeStyle, TextFace, faces, CircleFace
import argparse

args_nodes = ["population_size", "generation_time_in_year", "mutation_rate_per_generation"]


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
            radius = 15 * (float(getattr(node, arg)) - min_arg) / (max_arg - min_arg) + 5
        circle = CircleFace(radius=radius, color="RoyalBlue", style="sphere")
        circle.opacity = 0.3
        faces.add_face_to_node(circle, node, 0, position="float")
        for col, align in enumerate(args_nodes):
            faces.add_face_to_node(TextFace(getattr(node, align) + " "), node, col, position="aligned")


def tree_plot(input_simu):
    for arg in args_nodes:
        t = Tree(input_simu, format=3)
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.complete_branch_lines_when_necessary = False
        max_arg = max([float(getattr(n, arg)) for n in t.traverse()])
        min_arg = min([float(getattr(n, arg)) for n in t.traverse()])
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
