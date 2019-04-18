#!python3
from plot_module import *
import argparse
import os


def plot_simulation(input_simu, render):
    t = Tree(input_simu, format=3)

    nodes_file = input_simu.replace(".nhx", ".nodes.tsv")
    if os.path.isfile(nodes_file):
        nodes = pd.read_csv(nodes_file, sep='\t')
        ss_nodes = ["Branch_dN", "Branch_dS", "Branch_dNdS"]
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
                for key in [ss, ss + "_pred", ss + "_flow"]:
                    values = [float(getattr(n, key)) for n in t.traverse() if key in n.features]
                    if len(values) > 0 and len(values) == len(cols):
                        plt.scatter(list(cols), values, label="Exome " + key)

                if "theta" in ss:
                    plt.yscale("log")
                plt.ylabel(ss)
                plt.xticks(rotation='vertical')
                plt.legend()
                plt.savefig("{0}.tsv.{1}.png".format(input_simu.replace(".nhx", ""), ss), bbox_inches='tight')
                plt.clf()
                plt.close("all")

    args_nodes = set()
    for node in t.traverse():
        args_nodes = args_nodes.union(node.features)

    branch_dict = {}
    for arg in args_nodes:
        values = [float(getattr(n, arg)) for n in t.traverse() if
                  arg in n.features and convertible_to_float(getattr(n, arg))]
        if len(values) > 1 and render:
            ts = TreeStyle()
            ts.show_leaf_name = True
            ts.complete_branch_lines_when_necessary = False
            max_arg = max(values)
            min_arg = min(values)
            ts.layout_fn = lambda x: layout(x, arg, min_arg, max_arg)
            nameF = TextFace(arg, fsize=7)
            nameF.rotation = -90
            ts.aligned_header.add_face(nameF, column=0)
            ts.title.add_face(TextFace("{0} in simulation".format(arg), fsize=20), column=0)
            t.render("{0}.{1}.png".format(input_simu, arg), tree_style=ts)
        if len(values) > 1:
            if ("Branch" in arg) and (("dNd" in arg) or ("LogNe" in arg)):
                branch_dict[arg] = values

    der_pop_size = [(np.log10(float(n.population_size)) - float(n.Branch_LogNe)) for n in t.traverse() if
                    not n.is_root()]
    plot_correlation("{0}.correlation.png".format(input_simu), branch_dict, {}, der_pop_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--render', dest='render', action='store_true')
    parser.set_defaults(render=False)
    parser.add_argument('-t', '--tree', required=True, type=str, default='',
                        dest="t", metavar="<tree>", help="The tree to be drawn")
    args = parser.parse_args()
    plot_simulation(args.t, args.render)
