#!python3
import argparse
import os


def create_experiment(experiment):
    exp_path = os.getcwd() + '/Experiments/' + experiment

    os.makedirs(exp_path, exist_ok=True)
    os.system('cp config.yaml {0}'.format(exp_path))
    os.system('cp Snakefile {0}'.format(exp_path))
    run_file = exp_path + "/snakeslurm.sh"
    with open(run_file, 'w') as w:
        w.write("#!/usr/bin/env bash\n")
        run = 'snakemake -j 99 --cluster "sbatch -J {0} -p normal -N 1 ' \
              '-o {1}/slurm.%x.%j.out -e {1}/slurm.%x.%j.err '.format(experiment, exp_path)
        run += '--cpus-per-task={params.threads} --mem={params.mem} -t {params.time}"\n'
        w.write(run)
    os.system("chmod 755 " + run_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-e', '--experiment', required=False, type=str, default="", dest="experiment")
    args = parser.parse_args()
    if args.experiment == "":
        import uuid

        args.experiment = str(uuid.uuid4())
    create_experiment(args.experiment)
