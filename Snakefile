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
SIMULATION_PARAMS = '--newick ' + TREE
SIMULATION_PARAMS += ' --preferences ' + PREFERENCES
SIMULATION_PARAMS += ' --nuc_matrix ' + copy_params(EXPERIMENT, ROOT, config['SIMULATION']['NUC_MATRIX'])
SIMULATION_PARAMS += ' --precision_matrix ' + copy_params(EXPERIMENT, ROOT, config['SIMULATION']['PRECISION_MATRIX'])
SIMULATION_PARAMS += ' --mutation_rate_per_generation {0}'.format(config['SIMULATION']['MUTATION_RATE_PER_GENERATION'])
SIMULATION_PARAMS += ' --root_age {0}'.format(config['SIMULATION']['ROOT_AGE'])
SIMULATION_PARAMS += ' --generation_time {0}'.format(config['SIMULATION']['GENERATION_TIME'])
SIMULATION_PARAMS += ' --beta {0}'.format(config['SIMULATION']['BETA'])
SIMULATION_PARAMS += ' --exon_size {0}'.format(config['SIMULATION']['EXON_SIZE'])
SIMULATION_PARAMS += ' --seed {0}'.format(config['SIMULATION']['SEED'])
if config['SIMULATION']['BRANCH_WISE_CORRELATION']:
    SIMULATION_PARAMS += ' --branch_wise_correlation'
if config['SIMULATION']['FIXED_POP_SIZES']:
    SIMULATION_PARAMS += ' --fix_pop_size'
if config['SIMULATION']['FIXED_MUT_RATE']:
    SIMULATION_PARAMS += ' --fix_mut_rate'
if config['SIMULATION']['FIXED_GENERATION_TIME']:
    SIMULATION_PARAMS += ' --fix_gen_time'

SIMUPOLY_PARAMS = '--population_size {0}'.format(config['SimuPoly']['POP_SIZE'])
SIMUPOLY_PARAMS += ' --sample_size {0}'.format(config['SimuPoly']['SAMPLE_SIZE'])

SIMUDIV_PARAMS = '--nbr_grid_step {0}'.format(config['SimuDiv']['NBR_GRID_STEP'])

SIMULATION_SIMUMODE_PARAM = {"SimuDiv": SIMUDIV_PARAMS, "SimuPoly": SIMUPOLY_PARAMS}

PLOT_BURN_IN = config['PLOT']['BURN_IN']

# Parameters for the inference
assert (PLOT_BURN_IN < config['INFERENCE']['POINTS'])
INFERENCE_PARAMS = '-t ' + TREE
INFERENCE_PARAMS += ' -u {0}'.format(config['INFERENCE']['POINTS'])
INFERENCE_PARAMS += ' --precision {0}'.format(config['INFERENCE']['PRECISION'])
if config['INFERENCE']['CLAMP_PREFERENCES']:
    INFERENCE_PARAMS += ' --profiles ' + copy_params(EXPERIMENT, ROOT, PREFERENCES)
else:
    INFERENCE_PARAMS += ' --ncat {0}'.format(config['INFERENCE']['NCAT'])
if config['INFERENCE']['CLAMP_GENERATION_TIME']:
    INFERENCE_PARAMS += ' --clamp_gen_time'
if config['INFERENCE']['CLAMP_POP_SIZES']:
    INFERENCE_PARAMS += ' --clamp_pop_sizes'
if config['INFERENCE']['CLAMP_NUC_MATRIX']:
    INFERENCE_PARAMS += ' --clamp_nuc_matrix'
if config['INFERENCE']['CLAMP_CORR_MATRIX']:
    INFERENCE_PARAMS += ' --clamp_corr_matrix'
if config['INFERENCE']['GEODESIC']:
    INFERENCE_PARAMS += ' --geodesic'

INFERENCE_CHAINS = config['INFERENCE_REPLICATE']['CHAINS']
INFERENCE_POLYMORPHISM = config['INFERENCE_REPLICATE']['POLYMORPHISM']
INFERENCE_POLYMORPHISM_PARAM = {True: ' -p', False: ''}
INFERENCE_MODELS = config['INFERENCE_REPLICATE']['MODELS']
INFERENCE_SIMULATORS = config['INFERENCE_REPLICATE']['SIMULATORS']

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

localrules: all, all_prefs, all_trace, plot_simulation, plot_trace, plot_profiles, make_bayescode, bayescode_model, make_simuevol

rule all:
    input:
         EXPERIMENT + '/all_prefs',
         EXPERIMENT + '/all_trace'

rule make_bayescode:
    output: touch(EXPERIMENT + "/make_bayescode")
    input: dir=EXPERIMENT + '/bayescode.version'
    params: compile="&& make clean && make" if COMPILE else ""
    log: out=EXPERIMENT + '/bayescode.stdout', err=EXPERIMENT + '/bayescode.stderr'
    shell:
         'cd {ROOT}/bayescode {params.compile} 2> {log.err} 1> {log.out}'

rule bayescode_model:
    output:
          model=EXPERIMENT + '/{model}mutsel',
          read=EXPERIMENT + '/read{model}mutsel'
    input: EXPERIMENT + "/make_bayescode"
    shell:
         'cd {ROOT}/bayescode && cp _build/{wildcards.model}mutsel {EXPERIMENT} && cp _build/read{wildcards.model}mutsel {EXPERIMENT}'

rule make_simuevol:
    output:
          EXPERIMENT + '/SimuDiv',
          EXPERIMENT + '/SimuPoly'
    input: dir=EXPERIMENT + '/SimuEvol.version'
    params: compile="&& make clean && make" if COMPILE else ""
    log: out=EXPERIMENT + '/SimuEvol.stdout', err=EXPERIMENT + '/SimuEvol.stderr'
    shell: 'cd {ROOT}/SimuEvol {params.compile} 2> {log.err} 1> {log.out} && cp build/SimuDiv {EXPERIMENT} && cp build/SimuPoly {EXPERIMENT}'

rule run_simulation:
    output: touch(EXPERIMENT + '/{simumode}_exp')
    input:
         exec=EXPERIMENT + '/{simumode}',
         config_core=EXPERIMENT + '/config.SIMULATION',
         config_pan=EXPERIMENT + '/config.' + '{simumode}'
    params:
          time="0-23:00", mem=5000, threads=1,
          pan=lambda w: SIMULATION_SIMUMODE_PARAM[w.simumode]
    benchmark: EXPERIMENT + "/benchmarks.simulation.{simumode}.tsv"
    log: out=EXPERIMENT + '/{simumode}_exp.stdout', err=EXPERIMENT + '/{simumode}_exp.stderr'
    shell: '{input.exec} {SIMULATION_PARAMS} {params.pan} --output {output} 2> {log.err} 1> {log.out}'

rule plot_simulation:
    output: plot=touch(EXPERIMENT + '/simu_{simumode}_plot')
    input:
         src=ROOT + "/scripts/plot_simulation.py",
         simu=rules.run_simulation.output
    log: out=EXPERIMENT + '/{simumode}_plot.stdout', err=EXPERIMENT + '/{simumode}_plot.stderr'
    shell: 'python3 {input.src} --tree {input.simu}.nhx'

rule run_inference:
    output: touch(EXPERIMENT + '/{simumode}_{model}_{polymorphism}_{chain}_run')
    input:
         exec=rules.bayescode_model.output.model,
         simu=rules.run_simulation.output,
         param_infer=EXPERIMENT + '/config.INFERENCE'
    params:
          time="4-00:00", mem=5000, threads=1,
          poly=lambda w: INFERENCE_POLYMORPHISM_PARAM[w.polymorphism.lower() == 'true'],
          traits=lambda w, input: " --traitsfile {0}.traits.tsv".format(input.simu ) if (config['INFERENCE']['TRAITS'] and w.model == "node") else ""
    benchmark: EXPERIMENT + "/benchmarks.inference.{simumode}_{model}_{polymorphism}_{chain}_run.tsv"
    log: out=EXPERIMENT + '/{simumode}_{model}_{polymorphism}_{chain}_run.stdout', err=EXPERIMENT + '/{simumode}_{model}_{polymorphism}_{chain}_run.stderr'
    shell: '{input.exec} -a {input.simu}.ali {INFERENCE_PARAMS}{params.poly}{params.traits} {output} 2> {log.err} 1> {log.out}'

rule plot_trace:
    output: plot=directory(EXPERIMENT + '/inference_{simumode}_{model}_plot')
    input:
         src=ROOT + "/scripts/plot_trace.py",
         trace=expand(EXPERIMENT + '/{{simumode}}_{{model}}_{polymorphism}_{chain}_run', chain=INFERENCE_CHAINS, polymorphism=INFERENCE_POLYMORPHISM),
         simu=rules.run_simulation.output,
         param_plot=EXPERIMENT + '/config.PLOT'
    log: out=EXPERIMENT + '/{simumode}_{model}_inference_plot.stdout', err=EXPERIMENT + '/{simumode}_{model}_inference_plot.stderr'
    shell:
         'mkdir -p {output.plot} && python3 {input.src} --simu {input.simu} --trace {input.trace} --output {output.plot} --burn_in {PLOT_BURN_IN} 2> {log.err} 1> {log.out}'

rule all_trace:
    output: touch(EXPERIMENT + '/all_trace')
    input:
         expand(EXPERIMENT + '/simu_{simumode}_plot', simumode=INFERENCE_SIMULATORS),
         expand(EXPERIMENT + '/inference_{simumode}_{model}_plot', simumode=INFERENCE_SIMULATORS, model=INFERENCE_MODELS)

rule read_profiles:
    output: EXPERIMENT + '/{simumode}_{model}_{polymorphism}_{chain}_read.siteprofiles'
    input:
         trace=rules.run_inference.output,
         exec=rules.bayescode_model.output.read,
         param_plot=EXPERIMENT + '/config.PLOT'
    params: time="0-01:00", mem=5000, threads=1
    benchmark: EXPERIMENT + "/benchmarks.inference_{simumode}_{model}_{polymorphism}_{chain}_read.tsv"
    log: out=EXPERIMENT + '/{simumode}_{model}_{polymorphism}_{chain}_read.stdout', err=EXPERIMENT + '/{simumode}_{model}_{polymorphism}_{chain}_read.stderr'
    shell: '{input.exec} --burnin {PLOT_BURN_IN} -s --profiles {output} {input.trace} 2> {log.err} 1> {log.out}'

rule plot_profiles:
    output: prefs=directory(EXPERIMENT + '/{simumode}_{model}_inferred_prefs')
    input:
         src=ROOT + "/scripts/plot_profiles.py",
         profiles=expand(EXPERIMENT + '/{{simumode}}_{{model}}_{polymorphism}_{chain}_read.siteprofiles', chain=INFERENCE_CHAINS, polymorphism=INFERENCE_POLYMORPHISM)
    log: out=EXPERIMENT + '/{simumode}_{model}_inferred_prefs.stdout', err=EXPERIMENT + '/{simumode}_{model}_inferred_prefs.stderr'
    shell:
         'mkdir -p {output.prefs} && python3 {input.src} --input {PREFERENCES} --infer {input.profiles} --output {output.prefs} 2> {log.err} 1> {log.out}'

rule all_prefs:
    output: touch(EXPERIMENT + '/all_prefs')
    input:
         expand(EXPERIMENT + '/{simumode}_{model}_inferred_prefs', simumode=INFERENCE_SIMULATORS, model=INFERENCE_MODELS)
