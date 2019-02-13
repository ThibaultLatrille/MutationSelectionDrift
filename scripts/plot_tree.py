#!python3

from ete3 import Tree, TreeStyle, TextFace, faces, CircleFace
import argparse

args_nodes = ["population_size", "generation_time_in_year", "mutation_rate_per_generation"]


def layout(node, arg, min_arg, max_arg):
    if arg in node.features:
        radius = 15 * (float(getattr(node, arg)) - min_arg) / (max_arg - min_arg) + 15
        circle = CircleFace(radius=radius, color="RoyalBlue", style="sphere")
        circle.opacity = 0.3
        faces.add_face_to_node(circle, node, 0, position="float")
        for col, align in enumerate(args_nodes):
            faces.add_face_to_node(TextFace(getattr(node, align)), node, col, position="aligned")


def tree_plot(input_simu):
    for arg in args_nodes:
        t = Tree(input_simu)
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = True
        max_arg = max([float(getattr(n, arg)) for n in t.traverse()])
        min_arg = min([float(getattr(n, arg)) for n in t.traverse()])
        ts.layout_fn = lambda x: layout(x, arg, min_arg, max_arg)
        t.render("{0}.{1}.png".format(input_simu, arg), tree_style=ts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, default='',
                        dest="t", metavar="<tree>", help="The tree to be drawn")
    args = parser.parse_args()
    tree_plot(args.t)
