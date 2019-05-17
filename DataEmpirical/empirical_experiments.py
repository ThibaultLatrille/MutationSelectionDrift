#!python3
import argparse
import os
from subprocess import run


def create_experiment(name):
    for folder in os.listdir(os.getcwd()):
        if os.path.isdir(folder) and folder != "Experiments" and folder != "OrthoMam":
            experiment = name + "_" + folder
            exp_path = os.getcwd() + '/Experiments/' + experiment

            os.makedirs(exp_path, exist_ok=True)
            os.system('cp config.yaml {0}'.format(exp_path))
            os.system('cp ./{0}/CDS.ali {1}'.format(folder, exp_path))
            os.system('cp ./{0}/rootedtree.nhx {1}'.format(folder, exp_path))
            os.remove(exp_path + "/Snakefile") if os.path.exists(exp_path + "/Snakefile") else None
            os.symlink(os.getcwd() + "/Snakefile", exp_path + "/Snakefile")
            run_file = exp_path + "/snakeslurm.sh"
            with open(run_file, 'w') as w:
                w.write("#!/usr/bin/env bash\n")
                run_str = 'snakemake -j 99 --cluster "sbatch -J {0} -p normal -N 1 ' \
                          '-o {1}/slurm.%x.%j.out -e {1}/slurm.%x.%j.err '.format(experiment, exp_path)
                run_str += '--cpus-per-task={params.threads} --mem={params.mem} -t {params.time}"\n'
                w.write(run_str)
            os.system("chmod 755 " + run_file)
            screen = 'screen -dmS ' + experiment + ' bash -c "cd ' + exp_path + ' && ./snakeslurm.sh"'
            with open(exp_path + "/screen.sh", 'w') as w:
                w.write("#!/usr/bin/env bash\n")
                w.write(screen)
            print(screen)
            run(screen, shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-e', '--experiment', required=False, type=str, default="", dest="experiment")
    args = parser.parse_args()
    if args.experiment == "":
        import uuid

        args.experiment = str(uuid.uuid4())
    create_experiment(args.experiment)
