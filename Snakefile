import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from csv import reader
import os

def path_without_extension(path_to_file):
    return os.path.splitext(path_to_file)[0]

def cmd_to_stdout(cmd):
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

def print_c(txt):
    print('\033[94m' + txt + '\033[0m')


def copy_params(experiment, param_path):
    new_path = "{0}/{1}".format(experiment, os.path.basename(param_path))
    os.system('cp {0} {1}'.format(param_path, new_path))
    return new_path

def open_config(folder, config_file):
    config_path = '{0}/{1}'.format(folder, config_file)
    lock_path = path_without_extension(config_path) +  '.lock'
    if os.path.exists(config_path):
        diff = cmd_to_stdout('diff {0} {0}'.format(config_path, lock_path))
        if diff != '':
            print_c('Config file `{0}` has changed.'.format(config_path))
            print_c(diff)
            os.system('cp {0} {1}'.format(config_path, lock_path))
    else:
        os.system('cp {0} {1}'.format(config_file, config_path))
        os.system('cp {0} {1}'.format(config_path, lock_path))
    configfile: config_path



configfile: 'config.yaml'
EXPERIMENT = os.getcwd() + '/' + config['EXPERIMENT']
if os.path.exists(EXPERIMENT):
    print_c('Experiment {0} already exists.'.format(config['EXPERIMENT']))
else:
    os.mkdir(EXPERIMENT)
    print_c('Saving experiment into {0}'.format(EXPERIMENT))

os.system('cp Snakefile {0}/{1}.smk'.format(EXPERIMENT, config['EXPERIMENT']))

open_config(EXPERIMENT, 'params.yaml')

TREE = copy_params(EXPERIMENT, config['TREE'])
PREFERENCES = copy_params(EXPERIMENT, config['PREFERENCES'])

open_config(EXPERIMENT, 'params.simulation.yaml')
# Parameters for the simulation
SIMULATION = EXPERIMENT + '/simulation'
SIMULATION_PARAMS = '--newick ' + TREE
SIMULATION_PARAMS += ' --preferences ' + PREFERENCES
SIMULATION_PARAMS += ' --nuc_matrix ' + copy_params(EXPERIMENT, config['SIMULATION_NUC_MATRIX'])
SIMULATION_PARAMS += ' --correlation_matrix ' + copy_params(EXPERIMENT, config['SIMULATION_COR_MATRIX'])
SIMULATION_PARAMS += ' --mu {0}'.format(config['SIMULATION_MUTATION_RATE'])
SIMULATION_PARAMS += ' --root_age {0}'.format(config['SIMULATION_ROOT_AGE'])
SIMULATION_PARAMS += ' --generation_time {0}'.format(config['SIMULATION_GENERATION_TIME'])
SIMULATION_PARAMS += ' --pop_size {0}'.format(config['SIMULATION_POP_SIZE'])
SIMULATION_PARAMS += ' --sample_size {0}'.format(config['SIMULATION_SAMPLE_SIZE'])
SIMULATION_PARAMS += ' --beta {0}'.format(config['SIMULATION_BETA'])
if config['SIMULATION_LINKED_SITES']:
    SIMULATION_PARAMS += ' --linked'

open_config(EXPERIMENT, 'params.plot.yaml')
PLOT_BURN_IN = config['PLOT_BURN_IN']

open_config(EXPERIMENT, 'params.inference.yaml')
# Parameters for the inference
assert(PLOT_BURN_IN < config['INFERENCE_POINTS'])
INFERENCE = EXPERIMENT + '/inference'
INFERENCE_PARAMS = '-t ' + TREE
# INFERENCE_PARAMS += ' -c ' + copy_params(EXPERIMENT, config['INFERENCE_PREFERENCES'])
INFERENCE_PARAMS += ' --ncat {0}'.format(config['INFERENCE_NCAT'])
INFERENCE_PARAMS += ' -u {0}'.format(config['INFERENCE_POINTS'])

open_config(EXPERIMENT, 'config.inference.yaml')
INFERENCE_CHAINS = config['INFERENCE_CHAINS']
INFERENCE_POLYMORPHISM = config['INFERENCE_POLYMORPHISM']
INFERENCE_POLYMORPHISM_PARAM = {True: ' -p', False: ''}

if not os.path.exists('bayescode'):
    print_c('BayesCode directory not found, cloning the github repository.')
    os.system('git clone https://github.com/bayesiancook/bayescode')
if not os.path.exists('SimuEvol'):
    print_c('SimuEvol directory not found, cloning the github repository.')
    os.system('git clone https://github.com/ThibaultLatrille/SimuEvol')

for program in ['bayescode', 'SimuEvol']:
    stdout = cmd_to_stdout('cd {0} && git log -1'.format(program))
    stdout += cmd_to_stdout('cd {0} && git diff src'.format(program))
    version_file = EXPERIMENT + '/' + program + '.version'
    if (not os.path.exists(version_file)) or open(version_file, 'r').read() != stdout:
        if os.path.exists(version_file):
            print_c('{0} source code has changed, re-compiling and re-executing.'.format(program))
        with open(version_file, 'w') as version_append:
            version_append.write(stdout)

rule all:
    input:
        INFERENCE + '_plot'

rule make_bayescode:
    output:
        EXPERIMENT + '/datedmutsel'
    input:
        dir = EXPERIMENT + '/bayescode.version'
    log:
        out = EXPERIMENT + '/bayescode.stdout',
        err = EXPERIMENT + '/bayescode.stderr'
    shell:
        'cd bayescode && make clean && make 2> {log.err} 1> {log.out} && cp _build/datedmutsel {EXPERIMENT}'

rule make_simupoly:
    output:
        EXPERIMENT + '/SimuPoly'
    input:
        dir = EXPERIMENT + '/SimuEvol.version'
    log:
        out = EXPERIMENT + '/SimuEvol.stdout',
        err = EXPERIMENT + '/SimuEvol.stderr'
    shell:
        'cd SimuEvol && make clean && make 2> {log.err} 1> {log.out} && cp build/SimuPoly {EXPERIMENT}'

rule run_simulation:
    output:
        touch(SIMULATION)
    input:
        exec = rules.make_simupoly.output,
        param_simu = EXPERIMENT + '/params.simulation.lock'
    benchmark:
        EXPERIMENT + "/benchmarks.simulation.tsv"
    log:
        out = SIMULATION + '.stdout',
        err = SIMULATION + '.stderr'
    shell:
        '{input.exec} {SIMULATION_PARAMS} --output {output} 2> {log.err} 1> {log.out}'

rule run_inference:
    output:
        touch(INFERENCE + '_{polymorphism}_{chain}')
    input:
        exec = rules.make_bayescode.output,
        simu = rules.run_simulation.output,
        param_infer = EXPERIMENT + '/params.inference.lock'
    params:
        poly = lambda w: INFERENCE_POLYMORPHISM_PARAM[w.polymorphism.lower() == 'true']
    benchmark:
        EXPERIMENT + "/benchmarks.inference_{polymorphism}_{chain}.tsv"
    log:
        out = INFERENCE + '_{polymorphism}_{chain}.stdout',
        err = INFERENCE + '_{polymorphism}_{chain}.stderr'
    shell:
        '{input.exec} -a {input.simu}.ali {INFERENCE_PARAMS}{params.poly} {output} 2> {log.err} 1> {log.out}'

rule plot_trace:
    output:
        plot = touch(INFERENCE + '_plot')
    input:
        infer = expand(rules.run_inference.output, chain=INFERENCE_CHAINS, polymorphism=INFERENCE_POLYMORPHISM),
        param_plot = EXPERIMENT + '/params.plot.lock'
    run:
        traces = dict()
        for filename in input.infer:
            if os.path.exists(filename + '.trace'):
                with open(filename + '.trace', 'r') as trace_open:
                    trace = list(reader(trace_open, delimiter='\t'))
                    for i, param in enumerate(trace[0]):
                        if param not in traces:
                            traces[param] = dict()
                        traces[param][filename] = [float(line[i]) for line in trace[(1 + PLOT_BURN_IN):]]

        for param, traces_param in traces.items():
            my_dpi = 128
            fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
            for name, param_trace in sorted(traces_param.items(), key=lambda x: x[0]):
                style = "-"
                if "False" in name:
                    style = "--"
                plt.plot(range(len(param_trace)), param_trace, style, alpha=0.5, linewidth=1, label=name)
            plt.xlabel('Point')
            plt.ylabel(param)
            plt.legend()
            plt.tight_layout()
            plt.savefig('{0}_{1}.png'.format(output.plot, param), format='png')
            plt.clf()
            plt.close('all')
        open(output.plot, 'a').close()
