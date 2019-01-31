import os
from subprocess import run

cluster = False

if cluster:
    current_dir = "/panhome/tlatrill/SimuEvol"
    cmd = "qsub"
else:
    current_dir = "/home/thibault/SimuEvol"
    cmd = "sh"


def run_inference(protein):
    mut_bias = 1
    mu = 1e-6
    pop_size = 640
    beta = 2.5
    linked = "false"

    nbr_cpu = 1

    newick_path = "{0}/data_trees/{1}.newick".format(current_dir, protein)
    preferences_path = "{0}/data_prefs/{1}.txt".format(current_dir, protein)

    simu_cmd = current_dir + "/SimuPoly --preferences={0} --newick={1} --mu={2} --gen=0.1 --lambda={3} " \
                             "--pop_size={4} --sample_size={5} --beta={6} --linked={7} --output={8}\n"

    infer_cmd = current_dir + "/aamutsel -d {0} -T {1} -p -ncat 100 -x 1 400 {2}_{3}\n"

    for qsub_id, sample_size in enumerate([1, 2, 4, 8, 12, 16, 24, 32]):
        qsub_name = "{0}_{1}_{2}_{3}".format(sample_size, protein, qsub_id, pop_size)

        folder_path = "{0}/poly_sample_size_{1}".format(current_dir, qsub_name)
        os.makedirs(folder_path, exist_ok=True)

        qsub_path = "{0}/qsub/{1}.pbs".format(current_dir, qsub_name)

        qsub_str = "#!/bin/bash\n"
        qsub_str += "#\n"
        qsub_str += "#PBS -q q1day\n"
        qsub_str += "#PBS -l nodes=1:ppn={0},mem=4gb\n".format(nbr_cpu)
        qsub_str += "#PBS -o /pandata/tlatrill/read/read_out_{0}\n".format(qsub_name)
        qsub_str += "#PBS -e /pandata/tlatrill/read/read_err_{0}\n".format(qsub_name)
        qsub_str += "#PBS -j oe\n"
        qsub_str += "#PBS -W umask=022\n"
        qsub_str += "#PBS -r n\n"
        qsub_str += "#PBS -r n\n"

        # Run SimuPoly
        output_path = "{0}/{1}".format(folder_path, qsub_name)

        qsub_str += simu_cmd.format(preferences_path, newick_path, mu, mut_bias,
                                    pop_size, sample_size, beta, linked, output_path)

        qsub_str += infer_cmd.format(output_path + ".ali", newick_path, output_path, 1)

        qsub_str += infer_cmd.format(output_path + ".ali", newick_path, output_path, 2)

        # qsub_str += "rm -f {0}\n".format(qsub_path)

        qsub = open(qsub_path, 'w')
        qsub.write(qsub_str)
        qsub.close()

        print("Running " + qsub_path)
        # run("{0} {1}".format(cmd, qsub_path), shell=True)
        print("Finished running")


for p in ["gal4", "np"]:
    run_inference(p)
