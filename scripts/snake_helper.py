import matplotlib
from ete3 import Tree, TreeStyle, TextFace, faces, AttrFace, CircleFace

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from csv import reader
import os
import yaml

args_nodes = ["population_size", "generation_time_in_year", "mutation_rate_per_generation"]


class DictDiffer(object):
    """
    Calculate the difference between two dictionaries as:
    (1) items added
    (2) items removed
    (3) keys same in both but changed values
    (4) keys same in both and unchanged values
    """

    def __init__(self, current_dict, past_dict):
        self.current_dict, self.past_dict = current_dict, past_dict
        self.set_current, self.set_past = set(current_dict.keys()), set(past_dict.keys())
        self.intersect = self.set_current.intersection(self.set_past)

    def added(self):
        return self.set_current - self.intersect

    def removed(self):
        return self.set_past - self.intersect

    def changed(self):
        return set(o for o in self.intersect if self.past_dict[o] != self.current_dict[o])

    def unchanged(self):
        return set(o for o in self.intersect if self.past_dict[o] == self.current_dict[o])

    def diff(self):
        return len(self.added()) + len(self.removed()) + len(self.changed())

    def __repr__(self):
        str_ = ""
        for o in self.added():
            str_ += "    Added parameter {0}={1}\n".format(o, self.current_dict[o])
        for o in self.removed():
            str_ += "    Removed parameter {0}\n".format(o)
        for o in self.changed():
            str_ += "    Parameter {0}={1} (was {2} before)\n".format(o, self.current_dict[o], self.past_dict[o])
        return str_


def path_without_extension(path_to_file):
    return os.path.splitext(path_to_file)[0]


def cmd_to_stdout(cmd):
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')


def print_c(txt):
    print('\033[94m' + txt + '\033[0m')


def open_yaml(yaml_path):
    return yaml.load(open(yaml_path, 'r'))


def copy_params(experiment, param_path):
    new_path = "{0}/{1}".format(experiment, os.path.basename(param_path))
    os.system('cp {0} {1}'.format(param_path, new_path))
    return new_path


def open_config(folder, config_file):
    config_path = '{0}/{1}'.format(folder, config_file)
    lock_path = path_without_extension(config_path) + '.lock'
    if os.path.exists(config_path):
        config_dict = open_yaml(config_path)
        lock_dict = open_yaml(lock_path)
        for step, step_dict in config_dict.items():
            dict_diff = DictDiffer(step_dict, lock_dict[step])
            if dict_diff.diff() > 0:
                print_c(step + " config has changed:")
                print_c(str(dict_diff))
                os.system("touch {0}/config.{1}".format(folder, step))
                os.system('cp {0} {1}'.format(config_path, lock_path))
    else:
        os.system('cp {0} {1}'.format(config_file, config_path))
        os.system('cp {0} {1}'.format(config_path, lock_path))
        config_dict = open_yaml(config_path)
    for step in config_dict:
        conf = "{0}/config.{1}".format(folder, step)
        if not os.path.exists(conf):
            os.system("touch " + conf)
    return config_path


def open_tsv(filepath):
    if os.path.exists(filepath):
        with open(filepath, 'r') as tsv_open:
            tsv = list(reader(tsv_open, delimiter='\t'))
        return tsv
    else:
        return [[]]


def to_float(element):
    try:
        return float(element)
    except ValueError:
        return element


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
        t = Tree(input_simu + ".nhx")
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = True
        max_arg = max([float(getattr(n, arg)) for n in t.traverse()])
        min_arg = min([float(getattr(n, arg)) for n in t.traverse()])
        ts.layout_fn = lambda x: layout(x, arg, min_arg, max_arg)
        t.render("{0}.{1}.png".format(input_simu, arg), tree_style=ts)


tree_plot("/home/thibault/PolyMutSel/Experiments/Tree_output/simulation")


def trace_plot(input_simu, input_infer, output_plot, burn_in):
    simu_params, simu_nodes, simu_taxa, traces = dict(), dict(), dict(), dict()

    for simu in input_simu:
        params = open_tsv(simu + '.parameters.tsv')
        for i, param in enumerate(params[0]):
            simu_params[param] = to_float(params[1][i])

        nodes = open_tsv(simu + '.nodes.tsv')
        header = nodes[0]
        index = {v: k for k, v in enumerate(header)}["index"]
        for line in nodes[1:]:
            simu_nodes[line[index]] = dict()
            for i, value in enumerate(line):
                simu_nodes[line[index]][header[i]] = to_float(value)

        taxa = open_tsv(simu + '.tsv')
        header = taxa[0]
        index = {v: k for k, v in enumerate(header)}["taxon_name"]
        for line in nodes[1:]:
            simu_nodes[line[index]] = dict()
            for i, value in enumerate(line):
                simu_nodes[line[index]][header[i]] = to_float(value)

    for filename in input_infer:
        trace = open_tsv(filename + '.trace')
        for i, param in enumerate(trace[0]):
            if param not in traces:
                traces[param] = dict()
            traces[param][filename] = [float(line[i]) for line in trace[(1 + burn_in):]]

    for param, traces_param in traces.items():
        my_dpi = 128
        fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        for name, param_trace in sorted(traces_param.items(), key=lambda x: x[0]):
            style = "-"
            if "False" in name:
                style = "--"
            if ("BranchPopSize" in param) and (max(param_trace) > 10):
                print("Issue with " + name + ", " + param)
                y_values = [min(i, 10) for i in param_trace]
                plt.plot(range(len(param_trace)), y_values, style, alpha=0.5, linewidth=1, label="!" + name)
            else:
                plt.plot(range(len(param_trace)), param_trace, style, alpha=0.5, linewidth=1, label=name)

        if param in simu_params:
            # If it can be found in the simulation parameters
            plt.axhline(y=simu_params[param], xmin=0.0, xmax=1.0, color='r', label="Simulation")

        if param[0] == "*":
            #  Then it's a branch
            pass
        if param[0] == "@":
            #  Then it's a taxon
            pass
        plt.xlabel('Point')
        plt.ylabel(param)
        plt.legend()
        plt.tight_layout()
        plt.savefig('{0}_{1}.png'.format(output_plot, param), format='png')
        plt.clf()
        plt.close('all')
    open(output_plot, 'a').close()
