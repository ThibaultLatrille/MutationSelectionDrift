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
    return yaml.load(open(yaml_path, 'r'))


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

EXPERIMENT = os.path.abspath('.')
ROOT = os.path.abspath('../..')

configfile: open_config(EXPERIMENT, 'config.yaml')

COMPILE = config['EXPERIMENT']['COMPILE']
TREE = copy_params(EXPERIMENT, ROOT, config['SIMULATION']['TREE'])
PREFERENCES = copy_params(EXPERIMENT, ROOT, config['SIMULATION']['PREFERENCES'])

# Parameters for the simulation
SIMULATION = EXPERIMENT + '/simulation'
SIMULATION_PARAMS = '--newick ' + TREE
SIMULATION_PARAMS += ' --preferences ' + PREFERENCES
SIMULATION_PARAMS += ' --nuc_matrix ' + copy_params(EXPERIMENT, ROOT, config['SIMULATION']['NUC_MATRIX'])
SIMULATION_PARAMS += ' --correlation_matrix ' + copy_params(EXPERIMENT, ROOT, config['SIMULATION']['COR_MATRIX'])
SIMULATION_PARAMS += ' --mu {0}'.format(config['SIMULATION']['MUTATION_RATE'])
SIMULATION_PARAMS += ' --root_age {0}'.format(config['SIMULATION']['ROOT_AGE'])
SIMULATION_PARAMS += ' --generation_time {0}'.format(config['SIMULATION']['GENERATION_TIME'])
SIMULATION_PARAMS += ' --pop_size {0}'.format(config['SIMULATION']['POP_SIZE'])
SIMULATION_PARAMS += ' --sample_size {0}'.format(config['SIMULATION']['SAMPLE_SIZE'])
SIMULATION_PARAMS += ' --beta {0}'.format(config['SIMULATION']['BETA'])
if config['SIMULATION']['LINKED_SITES']:
    SIMULATION_PARAMS += ' --linked'

PLOT_BURN_IN = config['PLOT']['BURN_IN']

# Parameters for the inference
assert (PLOT_BURN_IN < config['INFERENCE']['POINTS'])
INFERENCE = EXPERIMENT + '/inference'
INFERENCE_PARAMS = '-t ' + TREE
if config['INFERENCE']['PREFERENCES_CLAMP']:
    INFERENCE_PARAMS += ' -c ' + copy_params(EXPERIMENT, ROOT, PREFERENCES)
INFERENCE_PARAMS += ' --ncat {0}'.format(config['INFERENCE']['NCAT'])
INFERENCE_PARAMS += ' --precision {0}'.format(config['INFERENCE']['PRECISION'])
INFERENCE_PARAMS += ' -u {0}'.format(config['INFERENCE']['POINTS'])

INFERENCE_CHAINS = config['INFERENCE_REPLICATE']['CHAINS']
INFERENCE_POLYMORPHISM = config['INFERENCE_REPLICATE']['POLYMORPHISM']
INFERENCE_POLYMORPHISM_PARAM = {True: ' -p', False: ''}

for program in ['bayescode', 'SimuEvol']:
    stdout = cmd_to_stdout('cd {0}/{1} && git log -1'.format(ROOT, program))
    stdout += cmd_to_stdout('cd {0}/{1} && git diff src utils'.format(ROOT, program))
    version_file = EXPERIMENT + '/' + program + '.version'
    if (not os.path.exists(version_file)) or open(version_file, 'r').read() != stdout:
        if os.path.exists(version_file):
            if COMPILE:
                print_c('{0} source code has changed, re-compiling and re-executing.'.format(program))
            else:
                print_c(
                    '{0} source code has changed, re-executing (but not recompiling, make sure it was recompiled before).'.format(
                        program))
        with open(version_file, 'w') as version_append:
            version_append.write(stdout)

localrules:all, plot_tree, plot_trace, plot_profiles, make_bayescode, make_simupoly

rule all:
    input: SIMULATION + '_plot', INFERENCE + '_plot', INFERENCE + '_prefs'

rule make_bayescode:
    output:
          dated=EXPERIMENT + '/datedmutsel',
          read=EXPERIMENT + '/readdatedmutsel'
    input: dir=EXPERIMENT + '/bayescode.version'
    params: compile=" make clean && make 2> {log.err} 1> {log.out} &&" if COMPILE else ""
    log: out=EXPERIMENT + '/bayescode.stdout', err=EXPERIMENT + '/bayescode.stderr'
    shell:
         'cd {ROOT}/bayescode &&{params.compile} cp _build/datedmutsel {EXPERIMENT} && cp _build/readdatedmutsel {EXPERIMENT}'

rule make_simupoly:
    output: EXPERIMENT + '/SimuPoly'
    input: dir=EXPERIMENT + '/SimuEvol.version'
    params: compile=" make clean && make 2> {log.err} 1> {log.out} &&" if COMPILE else ""
    log: out=EXPERIMENT + '/SimuEvol.stdout', err=EXPERIMENT + '/SimuEvol.stderr'
    shell: 'cd {ROOT}/SimuEvol &&{params.compile} cp build/SimuPoly {EXPERIMENT}'

rule run_simulation:
    output: touch(SIMULATION)
    input:
         exec=rules.make_simupoly.output,
         param_simu=EXPERIMENT + '/config.SIMULATION'
    params: time="0-23:00", mem=5000, threads=1
    benchmark: EXPERIMENT + "/benchmarks.simulation.tsv"
    log: out=SIMULATION + '.stdout', err=SIMULATION + '.stderr'
    shell: '{input.exec} {SIMULATION_PARAMS} --output {output} 2> {log.err} 1> {log.out}'

rule plot_tree:
    output: plot=touch(SIMULATION + '_plot')
    input:
         src=ROOT + "/scripts/plot_tree.py",
         simu=rules.run_simulation.output
    shell: 'python3 {input.src} --tree {input.simu}.nhx'

rule run_inference:
    output: touch(INFERENCE + '_{polymorphism}_{chain}_run')
    input:
         exec=rules.make_bayescode.output.dated,
         simu=rules.run_simulation.output,
         param_infer=EXPERIMENT + '/config.INFERENCE'
    params:
          time="2-00:00", mem=5000, threads=1,
          poly=lambda w: INFERENCE_POLYMORPHISM_PARAM[w.polymorphism.lower() == 'true']
    benchmark: EXPERIMENT + "/benchmarks.inference_{polymorphism}_{chain}_run.tsv"
    log: out=INFERENCE + '_{polymorphism}_{chain}_run.stdout', err=INFERENCE + '_{polymorphism}_{chain}_run.stderr'
    shell: '{input.exec} -a {input.simu}.ali {INFERENCE_PARAMS}{params.poly} {output} 2> {log.err} 1> {log.out}'

rule read_profiles:
    output: INFERENCE + '_{polymorphism}_{chain}_read.siteprofiles'
    input:
         trace=rules.run_inference.output,
         exec=rules.make_bayescode.output.read,
         param_plot=EXPERIMENT + '/config.PLOT'
    params: time="0-01:00", mem=5000, threads=1
    benchmark: EXPERIMENT + "/benchmarks.inference_{polymorphism}_{chain}_read.tsv"
    log: out=INFERENCE + '_{polymorphism}_{chain}_read.stdout', err=INFERENCE + '_{polymorphism}_{chain}_read.stderr'
    shell: '{input.exec} --burnin {PLOT_BURN_IN} -s --profiles {output} {input.trace} 2> {log.err} 1> {log.out}'

rule plot_profiles:
    output: prefs=directory(INFERENCE + '_prefs')
    input:
         src=ROOT + "/scripts/plot_profiles.py",
         profiles=expand(rules.read_profiles.output, chain=INFERENCE_CHAINS, polymorphism=INFERENCE_POLYMORPHISM)
    log: out=INFERENCE + '_prefs.stdout', err=INFERENCE + '_prefs.stderr'
    shell:
         'mkdir -p {output.prefs} && python3 {input.src} --input {PREFERENCES} --infer {input.profiles} --output {output.prefs} 2> {log.err} 1> {log.out}'

rule plot_trace:
    output: plot=directory(INFERENCE + '_plot')
    input:
         src=ROOT + "/scripts/plot_trace.py",
         trace=expand(rules.run_inference.output, chain=INFERENCE_CHAINS, polymorphism=INFERENCE_POLYMORPHISM),
         simu=rules.run_simulation.output,
         param_plot=EXPERIMENT + '/config.PLOT'
    log: out=INFERENCE + '_plot.stdout', err=INFERENCE + '_plot.stderr'
    shell:
         'mkdir -p {output.plot} && python3 {input.src} --simu {input.simu} --trace {input.trace} --output {output.plot} --burn_in {PLOT_BURN_IN} 2> {log.err} 1> {log.out}'
