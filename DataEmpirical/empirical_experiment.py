#!python3
import argparse
import os
from subprocess import run


def create_experiment(prefix, name, cds_name, tree_name, lht, calibs, screen, sbatch, nbr_cpu):
    root_path = os.getcwd() + "/" + name
    experiment = prefix + "_{0}_{1}_{2}".format(name, cds_name, tree_name)
    exp_path = os.getcwd() + '/Experiments/' + experiment

    os.makedirs(exp_path, exist_ok=True)
    os.system('cp config.yaml {0}'.format(exp_path))
    os.system('cp {0}/{1} {2}/CDS.ali'.format(root_path, cds_name, exp_path))
    os.system('cp {0}/{1} {2}/rootedtree.nhx'.format(root_path, tree_name, exp_path))
    os.remove(exp_path + "/Snakefile") if os.path.exists(exp_path + "/Snakefile") else None
    os.symlink(os.getcwd() + "/Snakefile", exp_path + "/Snakefile")

    if os.path.isfile('{0}/{1}'.format(root_path, lht)):
        print("Life-History-Traits file provided (" + lht + ")")
        os.system('cp {0}/{1} {2}/life_history_traits.tsv'.format(root_path, lht, exp_path))

    if os.path.isfile('{0}/{1}'.format(root_path, calibs)):
        print("Fossil Calibrations file provided (" + calibs + ")")
        os.system('cp {0}/{1} {2}/calibs.tsv'.format(root_path, calibs, exp_path))

    if os.path.isfile('{0}/known_population_size.tsv'.format(root_path)):
        print("Known population size file provided (known_population_size.tsv)")
        os.system('cp {0}/known_population_size.tsv {1}'.format(root_path, exp_path))

    run_file = exp_path + "/snakeslurm.sh"
    with open(run_file, 'w') as w:
        w.write("#!/usr/bin/env bash\n")
        run_str = 'snakemake '
        if sbatch:
            run_str += '-j 99 --cluster "sbatch -J {0} -p long -N 1 ' \
                       '-o {1}/slurm.%x.%j.out -e {1}/slurm.%x.%j.err '.format(experiment, exp_path)
            run_str += '--cpus-per-task={params.threads} --mem={params.mem} -t {params.time}"\n'
        else:
            run_str += "-j {0} --printshellcmds".format(nbr_cpu)
        w.write(run_str)
    os.system("chmod 755 " + run_file)
    cmd = 'cd ' + exp_path + ' && ./snakeslurm.sh'
    screen_cmd = 'screen -dmS ' + "{0}_{1}".format(prefix, name) + ' bash -c "' + cmd + '"'
    with open(exp_path + "/screen.sh", 'w') as w:
        w.write("#!/usr/bin/env bash\n")
        w.write(screen_cmd)
    if screen:
        print(screen_cmd)
        run(screen_cmd, shell=True)
    else:
        print(cmd)
        run(cmd, shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--prefix', required=False, type=str, default="Cat50", dest="prefix")
    parser.add_argument('-n', '--name', required=False, type=str, default="Primates", dest="name")
    parser.add_argument('--cds', required=False, type=str, default="CDS.ali", dest="cds")
    parser.add_argument('--tree', required=False, type=str, default="rootedtree.nhx", dest="tree")
    parser.add_argument('--lht', required=False, type=str, default="life_history_traits.tsv", dest="lht")
    parser.add_argument('--calibs', required=False, type=str, default="calibs.tsv", dest="calibs")
    parser.add_argument('-s', '--screen', required=False, type=bool, default=False, dest="screen")
    parser.add_argument('-b', '--sbatch', required=False, type=bool, default=False, dest="sbatch")
    parser.add_argument('-c', '--nbr_cpu', required=False, type=int, default=4, dest="nbr_cpu")
    args = parser.parse_args()
    create_experiment(args.prefix, args.name, args.cds, args.tree, args.lht, args.calibs, args.screen, args.sbatch, args.nbr_cpu)
