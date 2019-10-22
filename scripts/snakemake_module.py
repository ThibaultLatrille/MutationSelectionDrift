import subprocess
import os
import yaml


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
    return yaml.load(open(yaml_path, 'r'), Loader=yaml.FullLoader)


def copy_params(experiment, root, param_path):
    new_path = "{0}/{1}".format(experiment, os.path.basename(param_path))
    if not os.path.exists(new_path):
        os.system('cp {0}/{1} {2}'.format(root, param_path, new_path))
    return new_path


def open_config(folder, config_file):
    config_path = '{0}/{1}'.format(folder, config_file)
    lock_path = path_without_extension(config_path) + '.lock'
    if os.path.exists(lock_path):
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
        os.system('cp {0} {1}'.format(config_path, lock_path))
        config_dict = open_yaml(config_path)
    for step in config_dict:
        conf = "{0}/config.{1}".format(folder, step)
        if not os.path.exists(conf):
            os.system("touch " + conf)
    return config_path


def diff_git_dir(git_dir, version_file):
    stdout = cmd_to_stdout('cd {0} && git log -1'.format(git_dir))
    stdout += cmd_to_stdout('cd {0} && git diff src utils'.format(git_dir))
    version_file = version_file + '.version'
    if (not os.path.exists(version_file)) or open(version_file, 'r').read() != stdout:
        if os.path.exists(version_file):
            print_c('{0} source code has changed, re-compiling and re-executing.'.format(git_dir))
        with open(version_file, 'w') as version_append:
            version_append.write(stdout)
