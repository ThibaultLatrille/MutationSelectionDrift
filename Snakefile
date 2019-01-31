configfile: "config.yaml"

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from csv import reader
import os

ROOT_DIR = os.getcwd() + "/"
EXPERIMENT = ROOT_DIR + config["EXPERIMENT"]
print("Experiment is saved into " + EXPERIMENT)
TREE = ROOT_DIR + config["TREE"]

# Parameters for the simulation
SIMULATION = EXPERIMENT + "/simulation"
SIMULATION_PARAMS = "--newick " + TREE
SIMULATION_PARAMS += " --preferences " + ROOT_DIR + config["SIMULATION_PREFERENCES"]
SIMULATION_PARAMS += " --nuc_matrix " + ROOT_DIR + config["SIMULATION_NUC_MATRIX"]
SIMULATION_PARAMS += " --correlation_matrix " + ROOT_DIR + config["SIMULATION_COR_MATRIX"]
SIMULATION_PARAMS += " --mu {0}".format(config["SIMULATION_MUTATION_RATE"])
SIMULATION_PARAMS += " --root_age {0}".format(config["SIMULATION_ROOT_AGE"])
SIMULATION_PARAMS += " --generation_time {0}".format(config["SIMULATION_GENERATION_TIME"])
SIMULATION_PARAMS += " --pop_size {0}".format(config["SIMULATION_POP_SIZE"])
SIMULATION_PARAMS += " --sample_size {0}".format(config["SIMULATION_SAMPLE_SIZE"])
SIMULATION_PARAMS += " --beta {0}".format(config["SIMULATION_BETA"])
if config["SIMULATION_LINKED_SITES"]:
    SIMULATION_PARAMS += " --linked"

# Parameters for the inference
INFERENCE = EXPERIMENT + "/inference"
INFERENCE_PARAMS = "-t " + TREE
INFERENCE_PARAMS += " -u {0}".format(config["INFERENCE_POINTS"])
INFERENCE_PARAMS += " --ncat {0}".format(config["INFERENCE_NCAT"])
INFERENCE_CHAINS = config["INFERENCE_CHAINS"]
INFERENCE_POLYMORPHISM = config["INFERENCE_POLYMORPHISM"]
INFERENCE_POLYMORPHISM_PARAM = {True: " -p", False: ""}

if not os.path.exists("bayescode"):
    print("BayesCode directory not found, cloning the github repository.")
    os.system("git clone https://github.com/bayesiancook/bayescode")
if not os.path.exists("SimuEvol"):
    print("SimuEvol directory not found, cloning the github repository.")
    os.system("git clone https://github.com/ThibaultLatrille/SimuEvol")

def cmd_to_stdout(cmd):
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

for program in ["bayescode", "SimuEvol"]:
    stdout = cmd_to_stdout("cd {0} && git log --oneline -1".format(program))
    stdout += cmd_to_stdout("cd {0} && git diff src".format(program))
    version_file = EXPERIMENT + "/" + program + ".version"
    if os.path.exists(version_file) and open(version_file, 'r').read() == stdout:
        print("{0} version has not changed".format(program))
    else:
        print("Writing " + version_file)
        with open(version_file, 'w') as version_append:
            version_append.write(stdout)

rule all:
    input:
        INFERENCE + '_plot'

rule copy_config:
    output:
        EXPERIMENT + '/config.yaml'
    shell:
        'cp config.yaml ' + EXPERIMENT

rule make_bayescode:
    output:
        'bayescode/_build/datedmutsel'
    input:
        dir = EXPERIMENT + '/bayescode.version'
    shell:
        'cd bayescode && make clean && make'

rule make_simupoly:
    output:
        'SimuEvol/build/SimuPoly'
    input:
        dir = EXPERIMENT + '/SimuEvol.version'
    shell:
        'cd SimuEvol && make clean && make'

rule run_simulation:
    output:
        SIMULATION
    input:
        exec = 'SimuEvol/build/SimuPoly',
        config = EXPERIMENT + '/config.yaml'
    shell:
        '{input.exec} {SIMULATION_PARAMS} --output {output} && touch {output}'

rule run_inference:
    output:
        INFERENCE + '_{polymorphism}_{chain}'
    input:
        exec = 'bayescode/_build/datedmutsel',
        simu = SIMULATION
    params:
        poly = lambda w: INFERENCE_POLYMORPHISM_PARAM[w.polymorphism.lower() == "true"]
    shell:
        '{input.exec} -a {input.simu}.ali {INFERENCE_PARAMS}{params.poly} {output} && touch {output}'

rule plot_trace:
    input:
        expand(INFERENCE + '_{polymorphism}_{chain}', chain=INFERENCE_CHAINS, polymorphism=INFERENCE_POLYMORPHISM)
    output:
        plot = INFERENCE + '_plot'
    run:
        traces = dict()
        for file_name in input:
            with open(file_name + ".trace", 'r') as trace_open:
                trace = list(reader(trace_open, delimiter='\t'))
                for i, param in enumerate(trace[0]):
                    if param not in traces:
                        traces[param] = dict()
                    traces[param][file_name] = [float(line[i]) for line in trace[1:]]

        for param, traces_param in traces.items():
            my_dpi = 128
            fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
            for name, param_trace in traces_param.items():
                plt.plot(range(len(param_trace)), param_trace, alpha=0.5, linewidth=1, label=name)
            plt.xlabel('Point')
            plt.ylabel(param)
            plt.legend()
            plt.tight_layout()
            plt.savefig("{0}_{1}.png".format(output.plot, param), format="png")
        open(output.plot, 'a').close()
