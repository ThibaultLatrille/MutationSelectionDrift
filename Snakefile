from scripts.snake_helper import *

configfile: 'config.yaml'
EXPERIMENT = os.getcwd() + '/Experiments/' + config['EXPERIMENT']['NAME']
if os.path.exists(EXPERIMENT):
    print_c('Experiment folder {0} already exists.'.format(config['EXPERIMENT']['NAME']))
else:
    os.mkdir(EXPERIMENT)
    print_c('Saving experiment into {0}'.format(EXPERIMENT))

os.system('cp Snakefile {0}'.format(EXPERIMENT))

configfile: open_config(EXPERIMENT, 'config.yaml')

TREE = copy_params(EXPERIMENT, config['SIMULATION']['TREE'])
PREFERENCES = copy_params(EXPERIMENT, config['SIMULATION']['PREFERENCES'])

# Parameters for the simulation
SIMULATION = EXPERIMENT + '/simulation'
SIMULATION_PARAMS = '--newick ' + TREE
SIMULATION_PARAMS += ' --preferences ' + PREFERENCES
SIMULATION_PARAMS += ' --nuc_matrix ' + copy_params(EXPERIMENT, config['SIMULATION']['NUC_MATRIX'])
SIMULATION_PARAMS += ' --correlation_matrix ' + copy_params(EXPERIMENT, config['SIMULATION']['COR_MATRIX'])
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
assert(PLOT_BURN_IN < config['INFERENCE']['POINTS'])
INFERENCE = EXPERIMENT + '/inference'
INFERENCE_PARAMS = '-t ' + TREE
if config['INFERENCE']['PREFERENCES_CLAMP']:
    INFERENCE_PARAMS += ' -c ' + copy_params(EXPERIMENT, PREFERENCES)
INFERENCE_PARAMS += ' --ncat {0}'.format(config['INFERENCE']['NCAT'])
INFERENCE_PARAMS += ' --precision {0}'.format(config['INFERENCE']['PRECISION'])
INFERENCE_PARAMS += ' -u {0}'.format(config['INFERENCE']['POINTS'])

INFERENCE_CHAINS = config['INFERENCE_REPLICATE']['CHAINS']
INFERENCE_POLYMORPHISM = config['INFERENCE_REPLICATE']['POLYMORPHISM']
INFERENCE_POLYMORPHISM_PARAM = {True: ' -p', False: ''}

for program in ['bayescode', 'SimuEvol']:
    stdout = cmd_to_stdout('cd {0} && git log -1'.format(program))
    stdout += cmd_to_stdout('cd {0} && git diff src utils'.format(program))
    version_file = EXPERIMENT + '/' + program + '.version'
    if (not os.path.exists(version_file)) or open(version_file, 'r').read() != stdout:
        if os.path.exists(version_file):
           print_c('{0} source code has changed, re-compiling and re-executing.'.format(program))
        with open(version_file, 'w') as version_append:
           version_append.write(stdout)

localrules: all, plot_tree, plot_trace, plot_profiles, make_bayescode, make_simupoly

rule all:
    input: SIMULATION + '_plot', INFERENCE + '_plot', INFERENCE + '_prefs'

rule make_bayescode:
    output:
        dated = EXPERIMENT + '/datedmutsel',
        read = EXPERIMENT + '/readdatedmutsel'
    input: dir = EXPERIMENT + '/bayescode.version'
    log: out = EXPERIMENT + '/bayescode.stdout', err = EXPERIMENT + '/bayescode.stderr'
    shell: 'cd bayescode && make clean && make 2> {log.err} 1> {log.out} && cp _build/datedmutsel {EXPERIMENT} && cp _build/readdatedmutsel {EXPERIMENT}'

rule make_simupoly:
    output: EXPERIMENT + '/SimuPoly'
    input: dir = EXPERIMENT + '/SimuEvol.version'
    log: out = EXPERIMENT + '/SimuEvol.stdout', err = EXPERIMENT + '/SimuEvol.stderr'
    shell: 'cd SimuEvol && make clean && make 2> {log.err} 1> {log.out} && cp build/SimuPoly {EXPERIMENT}'

rule run_simulation:
    output: touch(SIMULATION)
    input:
        exec = rules.make_simupoly.output,
        param_simu = EXPERIMENT + '/config.SIMULATION'
    params: time="0-23:00", mem=5, threads=1
    benchmark: EXPERIMENT + "/benchmarks.simulation.tsv"
    log: out = SIMULATION + '.stdout', err = SIMULATION + '.stderr'
    shell: '{input.exec} {SIMULATION_PARAMS} --output {output} 2> {log.err} 1> {log.out}'

rule plot_tree:
    output: plot = touch(SIMULATION + '_plot')
    input: simu = rules.run_simulation.output
    shell: 'python3 scripts/plot_tree.py --tree {input.simu}.nhx'

rule run_inference:
    output: touch(INFERENCE + '_{polymorphism}_{chain}_run')
    input:
        exec = rules.make_bayescode.output.dated,
        simu = rules.run_simulation.output,
        param_infer = EXPERIMENT + '/config.INFERENCE'
    params:
        time="2-00:00", mem=5, threads=1,
        poly = lambda w: INFERENCE_POLYMORPHISM_PARAM[w.polymorphism.lower() == 'true']
    benchmark: EXPERIMENT + "/benchmarks.inference_{polymorphism}_{chain}_run.tsv"
    log: out = INFERENCE + '_{polymorphism}_{chain}_run.stdout', err = INFERENCE + '_{polymorphism}_{chain}_run.stderr'
    shell: '{input.exec} -a {input.simu}.ali {INFERENCE_PARAMS}{params.poly} {output} 2> {log.err} 1> {log.out}'

rule read_profiles:
    output: INFERENCE + '_{polymorphism}_{chain}_read.siteprofiles'
    input:
        trace = rules.run_inference.output,
        exec = rules.make_bayescode.output.read,
        param_plot = EXPERIMENT + '/config.PLOT'
    params: time="0-01:00", mem=5, threads=1
    benchmark: EXPERIMENT + "/benchmarks.inference_{polymorphism}_{chain}_read.tsv"
    log: out = INFERENCE + '_{polymorphism}_{chain}_read.stdout', err = INFERENCE + '_{polymorphism}_{chain}_read.stderr'
    shell: '{input.exec} --burnin {PLOT_BURN_IN} -s --profiles {output} {input.trace} 2> {log.err} 1> {log.out}'

rule plot_profiles:
    output: prefs = directory(INFERENCE + '_prefs')
    input:
        src = "scripts/plot_profiles.py",
        profiles = expand(rules.read_profiles.output, chain=INFERENCE_CHAINS, polymorphism=INFERENCE_POLYMORPHISM)
    log: out = INFERENCE + '_prefs.stdout', err = INFERENCE + '_prefs.stderr'
    shell: 'mkdir -p {output.prefs} && python3 {input.src} --input {PREFERENCES} --infer {input.profiles} --output {output.prefs} 2> {log.err} 1> {log.out}'

rule plot_trace:
    output: plot = directory(INFERENCE + '_plot')
    input:
        src = "scripts/plot_trace.py",
        trace = expand(rules.run_inference.output, chain=INFERENCE_CHAINS, polymorphism=INFERENCE_POLYMORPHISM),
        simu = rules.run_simulation.output,
        param_plot = EXPERIMENT + '/config.PLOT'
    log: out = INFERENCE + '_plot.stdout', err = INFERENCE + '_plot.stderr'
    shell: 'mkdir -p {output.plot} && python3 {input.src} --simu {input.simu} --trace {input.trace} --output {output.plot} --burn_in {PLOT_BURN_IN} 2> {log.err} 1> {log.out}'