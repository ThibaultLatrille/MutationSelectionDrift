#!python3
import argparse
import os
from subprocess import run
import pandas as pd
from ete3 import Tree


def import_ali(filepath):
    ali_dict = dict()
    with open(filepath, 'r') as ali_file:
        next(ali_file)
        for line in ali_file:
            if line != "\n":
                name, seq = line.replace("  ", " ").replace("\n", "").split(" ")
                ali_dict[name] = seq
    return ali_dict


def export_ali(filepath, ali_dict):
    seq_sizes = set([len(v) for v in ali_dict.values()])
    assert (len(seq_sizes) == 1)
    with open(filepath, 'w') as ali_file:
        ali_file.write("{0} {1}\n".format(len(ali_dict), seq_sizes.pop()))
        ali_file.write("\n".join([" ".join(id_seq) for id_seq in ali_dict.items()]))


def create_experiment(exons, replicate):
    genes = pd.read_csv("./filtered.list")
    tree = Tree("./rootedtree.nhx", format=1)
    print("{0} leaves for the rooted tree".format(len(tree)))

    for rep in range(replicate):
        experiment = "OrthMam_Exons{0}_Replicates{1}_Id{2}".format(exons, replicate, rep)
        exp_path = os.getcwd() + '/Experiments/' + experiment
        os.makedirs(exp_path, exist_ok=True)
        os.system('cp config.yaml {0}'.format(exp_path))
        os.remove(exp_path + "/Snakefile") if os.path.exists(exp_path + "/Snakefile") else None
        os.symlink(os.getcwd() + "/Snakefile", exp_path + "/Snakefile")

        os.system('cp rootedtree.nhx {0}'.format(exp_path))

        alignments = []
        taxa = set(tree.get_leaf_names())
        for selected in genes.sample(exons).values:
            alignments.append(import_ali("{0}/singlegene_alignments/{1}.ali".format(os.getcwd(), selected[0])))
            taxa = taxa.union(alignments[-1].keys())

        merge_alignment = {k: "" for k in taxa}
        for alignment in alignments:
            seq_len_set = set([len(s) for s in alignment.values()])
            assert(len(seq_len_set) == 1)
            size = seq_len_set.pop()
            for taxon in taxa:
                if taxon in alignment:
                    merge_alignment[taxon] += alignment[taxon]
                else:
                    merge_alignment[taxon] += "-" * size

        export_ali(exp_path + "/CDS.ali", merge_alignment)

        trimmed_tree = tree.copy()
        trimmed_tree.prune(taxa, preserve_branch_length=True)
        trimmed_tree.write(outfile="{0}/rootedtree.nhx".format(exp_path))

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
    parser.add_argument('-e', '--exons', required=False, type=int, default=10, dest="exons")
    parser.add_argument('-r', '--replicate', required=False, type=int, default=10, dest="replicate")
    args = parser.parse_args()
    create_experiment(args.exons, args.replicate)
