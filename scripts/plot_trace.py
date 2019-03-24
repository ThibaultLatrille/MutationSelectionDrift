#!python3

from plot_module import *
import argparse
import os


def plot_trace(input_simu, input_trace, output_plot, burn_in):
    traces = dict()

    simu_params = {k: v[0] for k, v in pd.read_csv(input_simu + '.parameters.tsv', sep='\t').items()}
    simupoly = "population_size" in simu_params

    filenames = ["Simulation", "dNdN0_predicted_SIMU", "dNdN0_sequence_wise_predicted_SIMU",
                 "dNdN0_count_based_SIMU", "dNdS_event_based_SIMU", "dNdS_count_based_SIMU"]
    if simupoly:
        filenames += ["Watterson_SIMU", "WattersonSynonymous_SIMU", "Pairwise_SIMU", "PairwiseSynonymous_SIMU",
                      "dNdN0_event_based_SIMU", "dNdS_events_sum_SIMU"]
    else:
        filenames += ["dNdN0_flow_based_SIMU"]

    for filepath in input_trace:
        filenames.append(os.path.basename(filepath))
        for param, vals in pd.read_csv(filepath + '.trace', sep='\t').items():
            if param not in traces:
                traces[param] = dict()
            traces[param][os.path.basename(filepath)] = vals[burn_in:]

    tree = Tree(input_simu + ".nhx", format=3)
    is_ultrametric(tree)
    for node in tree.traverse():
        if "mutation_rate" in node.features:
            rate = float(node.mutation_rate)
        else:
            rate = max(float(node.mutation_rate_per_generation) / float(node.generation_time_in_year), 1e-30)
        node.add_feature("LogNodeNe.Simulation", np.log10(float(node.population_size)))
        node.add_feature("LogNodeMutRate.Simulation", np.log10(rate))
        if not node.is_root():
            time = float(node.get_distance(node.up)) / simu_params["tree_max_distance_to_root_in_year"]
            node.add_feature("BranchTime.Simulation", float(time))
            if "Branch_mutation_rate" in node.features:
                branch_rate = float(node.Branch_mutation_rate)
            else:
                branch_rate = max(
                    float(node.Branch_mutation_rate_per_generation) / float(node.Branch_generation_time_in_year), 1e-30)
            node.add_feature("LogBranchNe.Simulation", float(node.Branch_LogNe))
            node.add_feature("LogBranchMutRate.Simulation", np.log10(rate))
            node.add_feature("LogBranchLength.Simulation", np.log10(time * branch_rate))
            node.add_feature("BranchdNdS.dNdN0_predicted_SIMU", float(node.Branch_dNdN0_predicted))
            node.add_feature("BranchdNdS.dNdN0_sequence_wise_predicted_SIMU",
                             float(node.Branch_dNdN0_sequence_wise_predicted))
            node.add_feature("BranchdNdS.dNdN0_count_based_SIMU", float(node.Branch_dNdN0_count_based))
            node.add_feature("BranchdNdS.dNdS_event_based_SIMU", float(node.Branch_dNdS_event_based))
            node.add_feature("BranchdNdS.dNdS_count_based_SIMU", float(node.Branch_dNdS_count_based))
            if simupoly:
                node.add_feature("BranchdNdS.dNdS_events_sum_SIMU", float(node.Branch_dNdS))
                node.add_feature("BranchdNdS.dNdN0_event_based_SIMU", float(node.Branch_dNdN0_event_based))
            else:
                node.add_feature("BranchdNdS.dNdN0_flow_based_SIMU", float(node.Branch_dNdN0_flow_based))
        if node.is_leaf() and simupoly:
            node.add_feature("LeafTheta.Simulation",
                             4 * float(node.population_size) * float(node.mutation_rate_per_generation))
            node.add_feature("LeafTheta.Watterson_SIMU", float(node.theta_watterson))
            node.add_feature("LeafTheta.WattersonSynonymous_SIMU", float(node.theta_watterson_syn))
            node.add_feature("LeafTheta.Pairwise_SIMU", float(node.theta_pairwise))
            node.add_feature("LeafTheta.PairwiseSynonymous_SIMU", float(node.theta_pairwise_syn))

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
                    values = np.log10(np.exp(param_trace))
                    node.add_feature("LogNodeNe." + filename, np.mean(values))
                    node.add_feature("LogNodeNe_var." + filename, np.var(values))
                elif "*BranchPopSize" in param:
                    values = np.log10(param_trace)
                    node.add_feature("LogBranchNe." + filename, np.mean(values))
                    node.add_feature("LogBranchNe_var." + filename, np.var(values))
                elif "*NodeRate" in param:
                    values = np.log10(np.exp(param_trace))
                    node.add_feature("LogNodeMutRate." + filename, np.mean(values))
                    node.add_feature("LogNodeMutRate_var." + filename, np.var(values))
                elif "*BranchRate" in param:
                    values = np.log10(param_trace)
                    node.add_feature("LogBranchMutRate." + filename, np.mean(values))
                    node.add_feature("LogBranchMutRate_var." + filename, np.var(values))
                elif "*BranchTime" in param:
                    assert (not node.is_root())
                    node.add_feature("BranchTime." + filename, np.mean(param_trace))
                    node.add_feature("BranchTime_var." + filename, np.var(param_trace))
                elif "*BranchLength" in param:
                    assert (not node.is_root())
                    node.add_feature("LogBranchLength." + filename, np.mean(np.log10(param_trace)))
                    node.add_feature("LogBranchLength_var." + filename, np.var(np.log10(param_trace)))
                elif "*BranchdNdS" in param:
                    assert (not node.is_root())
                    node.add_feature("BranchdNdS." + filename, np.mean(param_trace))
                    node.add_feature("BranchdNdS_var." + filename, np.var(param_trace))
                elif "*Theta" in param:
                    assert (node.is_leaf())
                    node.add_feature("LeafTheta." + filename, np.mean(param_trace))
                    node.add_feature("LeafTheta_var." + filename, np.var(param_trace))
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

    for arg in ["BranchTime", "LeafTheta", "BranchdNdS", "LogBranchLength", "LogNodeMutRate", "LogBranchMutRate",
                "LogNodeNe", "LogBranchNe"]:
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

        if len(axis_dict) > 1:
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
