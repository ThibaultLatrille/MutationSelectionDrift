#!python3
import argparse
import os
from subprocess import run
import pandas as pd
from ete3 import Tree


def import_ali(filepath):
    ali_dict = dict()
    with open(filepath, 'r') as ali_file:
        sp, size = ali_file.readline().split(" ")
        for line in ali_file:
            if line != "\n":
                name, seq = line.replace("\t", " ").replace("  ", " ").replace("\n", "").split(" ")
                assert(int(size) == len(seq))
                ali_dict[name] = seq
    assert(int(sp) == len(ali_dict))
    return ali_dict


def export_ali(filepath, ali_dict):
    seq_sizes = set([len(v) for v in ali_dict.values()])
    assert (len(seq_sizes) == 1)
    with open(filepath, 'w') as ali_file:
        ali_file.write("{0} {1}\n".format(len(ali_dict), seq_sizes.pop()))
        ali_file.write("\n".join([" ".join(id_seq) for id_seq in ali_dict.items()]))
        for name, seq in ali_dict.items():
            for i, n in enumerate(seq):
                if n != "A" and n != "C" and n != "G" and n != "T" and n != "-" and n != "N" and n != "?":
                    print("Unexpected character {0} in sequence {1} at position {2}".format(n, name, i + 1))


def create_experiment(prefix, name, sample, replicate, tree_name, cds_list, lht, calibs, intersection, screen, sbatch, nbr_cpu, random_state):
    root_path = os.getcwd() + "/" + name
    tree = Tree("{0}/{1}".format(root_path, tree_name), format=1)
    print("{0} extant species found for the rooted tree/".format(len(tree)))

    if os.path.isfile('{0}/{1}'.format(root_path, cds_list)):
        print("Found list of CDS : " + cds_list)
        genes = pd.read_csv("{0}/{1}".format(root_path, cds_list), header=None)
    else:
        genes = pd.DataFrame(
            [i.replace(".ali", "") for i in os.listdir("{0}/singlegene_alignments".format(root_path)) if ".ali" in i])

    print("{0} CDS provided.".format(len(genes)))

    if replicate == -1:
        replicate = len(genes)

    for rep in range(replicate):
        random_state += 654
        experiment = prefix + "_{0}_{1}_{2}_Sample{3}_Replicates{4}_Id{5}".format(name, tree_name, cds_list, sample,
                                                                                 replicate, rep)
        exp_path = os.getcwd() + '/Experiments/' + experiment
        os.makedirs(exp_path, exist_ok=True)
        os.system('cp config.yaml {0}'.format(exp_path))
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

        alignments = []
        taxa = set(tree.get_leaf_names()) if intersection else set()

        if sample == -1:
            vals = [genes.loc[rep, :]]
        else:
            vals = genes.sample(sample, random_state=random_state).values
        for selected in vals:
            alignments.append(import_ali("{0}/singlegene_alignments/{1}.ali".format(root_path, selected[0])))
            taxa = taxa.intersection(alignments[-1].keys()) if intersection else taxa.union(alignments[-1].keys())

        taxa = taxa.intersection(set(tree.get_leaf_names()))

        merge_alignment = {k: "" for k in taxa}
        for alignment in alignments:
            seq_len_set = set([len(s) for s in alignment.values()])
            assert (len(seq_len_set) == 1)
            size = seq_len_set.pop()
            for taxon in taxa:
                if taxon in alignment:
                    merge_alignment[taxon] += alignment[taxon]
                else:
                    merge_alignment[taxon] += "-" * size

        export_ali(exp_path + "/CDS.ali", merge_alignment)

        trimmed_tree = tree.copy()
        trimmed_tree.prune(taxa, preserve_branch_length=True)
        print("{0} taxa for replicate {1}".format(len(taxa), rep + 1))
        trimmed_tree.write(outfile="{0}/rootedtree.nhx".format(exp_path), format=1)

        assert (set(merge_alignment.keys()) == set(
            Tree("{0}/rootedtree.nhx".format(exp_path), format=1).get_leaf_names()))

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
        screen_cmd = 'screen -dmS ' + "{0}_{1}_{2}".format(prefix, name, rep) + ' bash -c "' + cmd + '"'
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
    parser.add_argument('-n', '--name', required=True, type=str, default="OrthoMam", dest="name")
    # name can be one of ["Vertebrates", "47SP", "OrthoMam", "Isopods", "Primates"]
    parser.add_argument('--tree', required=True, type=str, default="rootedtree.lht.nhx", dest="tree")
    parser.add_argument('--cds', required=True, type=str, default="cds.highcoverage.list", dest="cds")
    parser.add_argument('--sample', required=False, type=int, default=18, dest="sample")
    parser.add_argument('--replicate', required=False, type=int, default=4, dest="replicate")
    parser.add_argument('--lht', required=False, type=str, default="", dest="lht")
    parser.add_argument('--calibs', required=False, type=str, default="", dest="calibs")
    parser.add_argument('--intersection', required=False, type=bool, default=False, dest="intersection")
    parser.add_argument('-s', '--screen', required=False, type=bool, default=False, dest="screen")
    parser.add_argument('-b', '--sbatch', required=False, type=bool, default=False, dest="sbatch")
    parser.add_argument('-c', '--nbr_cpu', required=False, type=int, default=4, dest="nbr_cpu")
    parser.add_argument('--random_state', required=False, type=int, default=0, dest="random_state")

    args = parser.parse_args()
    create_experiment(args.prefix, args.name, args.sample, args.replicate, args.tree, args.cds, args.lht, args.calibs,
                      args.intersection, args.screen, args.sbatch, args.nbr_cpu, args.random_state)
