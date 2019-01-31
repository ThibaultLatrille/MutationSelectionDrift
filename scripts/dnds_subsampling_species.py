from ete3 import Tree
from Bio import SeqIO
import os
from subprocess import run

cluster = False
nbr_cpu = 8
burn_in = 1
chain_length = 15
chains = ["1"]
min_n = 25

if cluster:
    current_dir = "/panhome/tlatrill/SimuEvol"
    pb_mpi_path = "/panhome/tlatrill/pbmpi2/data"

else:
    current_dir = "/home/thibault/SimuEvol"
    pb_mpi_path = "/home/thibault/Tools/pbmpi2/data"

protein = "np"

# Run SimuEvol for the 3 conditions
newick_path = "{0}/data_trees/{1}.newick".format(current_dir, protein)
prefs_path = "{0}/data_prefs/{1}.txt".format(current_dir, protein)
data_path = "{0}/data_pb".format(current_dir)

os.makedirs(current_dir + "/qsub", exist_ok=True)
os.makedirs(data_path, exist_ok=True)

for s, p in [(0.0, 1.0), (10.0, 0.0), (0.0, 0.0)]:
    output_path = "{0}/data_alignment/{1}_{2}_{3}".format(current_dir, protein, s, p)

    simu_evol_cmd = current_dir + "/SimuEvol --preferences={0} --newick={1} --output={2} --mu={3} --lambda={4}"
    simu_evol_cmd = simu_evol_cmd.format(prefs_path, newick_path, output_path, 2.5, 1.0)
    simu_evol_cmd += " --s={0} --p={1} --a=False".format(s, p)

    print("Running SimuEvol")
    print(simu_evol_cmd)
    run(simu_evol_cmd, shell=True)
    print("Finished running SimuEvol")

    fasta_seqs = SeqIO.parse(open(output_path + ".fasta", 'r'), 'fasta')
    ids_seqs = [(fasta.id, str(fasta.seq)) for fasta in fasta_seqs]

    len_set = set([len(seq) for id_seq, seq in ids_seqs])
    assert len(len_set) == 1, "Sequences have different length"
    len_seq = len_set.pop()
    assert len_seq % 3 == 0, "Sequences are not multiple of 3"

    t = Tree(newick_path)

    for index, node in enumerate(t.traverse("levelorder")):
        leaves = [l.name for l in node.get_leaves()]
        n = len(leaves)
        if n > min_n:
            ali_path = "{0}/data_alignment/{1}_{2}.ali".format(current_dir, protein, index)
            filtered_seq = [(id_seq, seq) for id_seq, seq in ids_seqs if (id_seq in leaves)]
            ali_file = open(ali_path, 'w')
            ali_file.write("{0} {1}\n".format(len(filtered_seq), len(filtered_seq[0][1])))
            ali_file.write("\n".join([" ".join(id_seq) for id_seq in filtered_seq]))
            ali_file.close()

            tree_path = "{0}/data_trees/{1}_{2}.newick".format(current_dir, protein, index)
            node.unroot()
            node.write(format=0, outfile=tree_path)

            # Create and submit .pbs for the subsampled trees
            for option, opt_name in [("-globalomega", "globalomega"),
                                     ("-siteomega", "siteomega"),
                                     ("-mutsel -dp", "mutsel"),
                                     ("-mutsel -freeomega -dp", "mutselfreeomega")]:
                for chain in chains:
                    run_id = "_".join(map(str, [protein, s, p, index, n, opt_name, chain]))

                    qsub_path = "{0}/qsub/{1}.pbs".format(current_dir, run_id)

                    qsub_str = "#!/bin/bash\n"
                    qsub_str += "#\n"
                    qsub_str += "#PBS -q q1day\n"
                    qsub_str += "#PBS -l nodes=1:ppn={0},mem=4gb\n".format(nbr_cpu)
                    qsub_str += "#PBS -o /pandata/tlatrill/out_err/out_{0}\n".format(run_id)
                    qsub_str += "#PBS -e /pandata/tlatrill/out_err/err_{0}\n".format(run_id)
                    qsub_str += "#PBS -j oe\n"
                    qsub_str += "#PBS -W umask=022\n"
                    qsub_str += "#PBS -r n\n"
                    qsub_str += "TMP=/tmp/tlatrill$RANDOM\n"
                    qsub_str += "export TMPDIR=$TMP\n"

                    qsub_str += "mpirun -n {0} {1}/pb_mpi -f -s -x 1 {2} {3}".format(nbr_cpu, pb_mpi_path, chain_length, option)
                    qsub_str += " -d {0}".format(ali_path)
                    qsub_str += " -T {0}".format(tree_path)
                    qsub_str += " {0}/{1}\n".format(data_path, run_id)

                    if "mutsel" in opt_name:
                        qsub_str += "{0}/readpb_mpi -x {1} -ss {2}/{3}\n".format(pb_mpi_path, burn_in, data_path, run_id)
                        qsub_str += "{0}/readpb_mpi -x {1} -om {2}/{3}\n".format(pb_mpi_path, burn_in, data_path, run_id)

                    qsub_str += "rm -rf $TMP\n"
                    qsub_str += "rm {0}\n".format(qsub_path)

                    qsub = open(qsub_path, 'w')
                    qsub.write(qsub_str)
                    qsub.close()

                    if cluster:
                        print("Submitting " + qsub_path + " to the cluster")
                        run("qsub {0}".format(qsub_path), shell=True)
                    else:
                        print("Running " + qsub_path)
                        run("sh {0}".format(qsub_path), shell=True)

print('Job completed')
