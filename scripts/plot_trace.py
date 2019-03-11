#!python3

from plot_module import *
import argparse
import os


def plot_trace(input_simu, input_trace, output_plot, burn_in):
    traces = dict()

    simu_params = {k: v[0] for k, v in pd.read_csv(input_simu + '.parameters.tsv', sep='\t').items()}

    filenames = ["Simulation", "Watterson_Simulation", "WattersonSynonymous_Simulation", "Pairwise_Simulation",
                 "PairwiseSynonymous_Simulation"]
    for filepath in input_trace:
        filenames.append(os.path.basename(filepath))
        for param, vals in pd.read_csv(filepath + '.trace', sep='\t').items():
            if param not in traces:
                traces[param] = dict()
            traces[param][os.path.basename(filepath)] = vals[burn_in:]

    tree = Tree(input_simu + ".nhx", format=3)
    is_ultrametric(tree)
    for node in tree.traverse():
        rate = max(float(node.mutation_rate_per_generation) / float(node.generation_time_in_year), 1e-30)
        node.add_feature("LogNe.Simulation", np.log(float(node.population_size)))
        node.add_feature("LogMutRate.Simulation", np.log(rate))
        if not node.is_root():
            time = float(node.get_distance(node.up)) / simu_params["tree_max_distance_to_root_in_year"]
            node.add_feature("Time.Simulation", float(time))
            node.add_feature("LogLength.Simulation", np.log(time * rate))
            node.add_feature("dNdS.Simulation", float(node.dNdS))
        if node.is_leaf():
            node.add_feature("Theta.Simulation",
                             4 * float(node.population_size) * float(node.mutation_rate_per_generation))
            node.add_feature("Theta.Watterson_Simulation", float(node.theta_watterson))
            node.add_feature("Theta.WattersonSynonymous_Simulation", float(node.theta_watterson_syn))
            node.add_feature("Theta.Pairwise_Simulation", float(node.theta_pairwise))
            node.add_feature("Theta.PairwiseSynonymous_Simulation", float(node.theta_pairwise_syn))
    for param, traces_param in traces.items():
        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        for filename, param_trace in sorted(traces_param.items(), key=lambda x: x[0]):
            style = "-"
            if "False" in filename:
                style = "--"
            plt.plot(range(len(param_trace)), param_trace, style, alpha=0.5, linewidth=1, label=filename)

            if param[0] == "*":
                node_name = param.split("_")[-1]

                if node_name != "":
                    node = next(tree.iter_search_nodes(name=node_name))
                else:
                    node = tree.get_tree_root()

                if "*NodePopSize" in param:
                    node.add_feature("LogNe." + filename, np.mean(param_trace))
                    node.add_feature("LogNe_var." + filename, np.var(param_trace))
                elif "*NodeRate" in param:
                    node.add_feature("LogMutRate." + filename, np.mean(param_trace))
                    node.add_feature("LogMutRate_var." + filename, np.var(param_trace))
                elif "*BranchTime" in param:
                    assert (not node.is_root())
                    node.add_feature("Time." + filename, np.mean(param_trace))
                    node.add_feature("Time_var." + filename, np.var(param_trace))
                elif "*BranchLength" in param:
                    assert (not node.is_root())
                    node.add_feature("LogLength." + filename, np.mean(np.log(param_trace)))
                    node.add_feature("LogLength_var." + filename, np.var(np.log(param_trace)))
                elif "*BranchdNdS" in param:
                    assert (not node.is_root())
                    node.add_feature("dNdS." + filename, np.mean(param_trace))
                    node.add_feature("dNdS_var." + filename, np.var(param_trace))
                elif "*Theta" in param:
                    assert (node.is_leaf())
                    node.add_feature("Theta." + filename, np.mean(param_trace))
                    node.add_feature("Theta_var." + filename, np.var(param_trace))
        if param.lower() in [s.lower() for s in simu_params]:
            # If it can be found in the Simulation parameters
            plt.axhline(y=simu_params[param], xmin=0.0, xmax=1.0, color='r', label="Simulation")

        plt.xlabel('Point')
        plt.ylabel(param)
        plt.legend()
        plt.tight_layout()
        plt.savefig('{0}/trace.{1}.png'.format(output_plot, param), format='png')
        plt.clf()
        plt.close('all')

    for arg in ["Time", "Theta", "dNdS", "LogLength", "LogMutRate", "LogNe"]:
        axis_dict = dict()
        err_dict = dict()
        for filename in filenames:
            attr = arg + "." + filename
            var = arg + "_var." + filename
            axis = np.array([getattr(node, attr) for node in tree.traverse() if attr in node.features])
            err = np.array([getattr(node, var) for node in tree.traverse() if var in node.features])
            if len(axis) > 0:
                axis_dict[filename] = axis
                if len(err) > 0:
                    err_dict[filename] = np.sqrt(err)

        for filename, axis in axis_dict.items():
            max_arg, min_arg = max(axis), min(axis)
            ts = TreeStyle()
            ts.show_leaf_name = True
            ts.complete_branch_lines_when_necessary = False
            columns = sorted(axis_dict)
            ts.layout_fn = lambda x: mutiple_layout(x, arg, min_arg, max_arg, filename, columns)
            for col, name in enumerate(columns):
                nameF = TextFace(name, fsize=7)
                nameF.rotation = -90
                ts.aligned_header.add_face(nameF, column=col)
            ts.title.add_face(TextFace("{1} in {2}".format(output_plot, arg, filename), fsize=20), column=0)
            tree.render("{0}/nhx.{1}.{2}.png".format(output_plot, arg, filename), tree_style=ts)

        plot_correlation('{0}/correlation.{1}.png'.format(output_plot, arg), axis_dict, err_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--simu', required=True, type=str, dest="simu")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-t', '--trace', required=True, type=str, nargs='+', dest="trace")
    parser.add_argument('-b', '--burn_in', required=False, type=int, default=0, dest="burn_in")
    args = parser.parse_args()
    plot_trace(args.simu, args.trace, args.output, args.burn_in)
