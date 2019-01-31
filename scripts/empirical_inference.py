import os
from subprocess import run

cluster = False

if cluster:
    current_dir = "/panhome/tlatrill/SimuEvol"
    cmd = "qsub"
else:
    current_dir = "/home/thibault/SimuEvol"
    cmd = "sh"

# folder = "data_empirical"
folder = "data_primates"
folder_path = "{0}/{1}".format(current_dir, folder)
fasta_files = [f for f in os.listdir(folder_path) if ".fasta" == f.strip()[-6:]]
print("Found {0} fasta files".format(len(fasta_files)))

for fasta_file in sorted(fasta_files):
    print("Run for fasta file {0}".format(fasta_file))
    prefix = fasta_file.strip().replace(".fasta", "")
    tree_path = "{0}/{1}.newick".format(folder_path, prefix)
    if not os.path.isfile(tree_path):
        print("Can't find the newick file {0}".format(tree_path))
    else:
        nbr_cpu = 4

        qsub_path = "{0}/{1}.pbs".format(folder_path, prefix)

        qsub_str = "#!/bin/bash\n"
        qsub_str += "#\n"
        qsub_str += "#PBS -q q1day\n"
        qsub_str += "#PBS -l nodes=1:ppn={0},mem=4gb\n".format(nbr_cpu)
        qsub_str += "#PBS -o /pandata/tlatrill/read/read_out_{0}\n".format(prefix)
        qsub_str += "#PBS -e /pandata/tlatrill/read/read_err_{0}\n".format(prefix)
        qsub_str += "#PBS -j oe\n"
        qsub_str += "#PBS -W umask=022\n"
        qsub_str += "#PBS -r n\n"
        qsub_str += "#PBS -r n\n"

        fasta_path = "{0}/{1}".format(folder_path, fasta_file)
        scripts_dir = "{0}/scripts".format(current_dir)

        rate = 0  # can be 0, 1, or 5
        freq = 3  # can be 1 or 3
        for omega in [1, 4, 95]:
            param = "{0}-{1}-{2}".format(rate, freq, omega)
            hyphy_batch_path = "{0}/{1}_{2}.bf".format(folder_path, prefix, param)
            batchfile_cmd = "python3 {0}/projected_mut_sel.py -d {0} -b {1} -f {2} -t {3} -p '{4}' \n"
            qsub_str += batchfile_cmd.format(scripts_dir, hyphy_batch_path, fasta_path, tree_path, param)
            qsub_str += "HYPHYMP {0} CPU={1}\n".format(hyphy_batch_path, nbr_cpu)

        qsub_str += "rm -f {0}\n".format(qsub_path)

        qsub = open(qsub_path, 'w')
        qsub.write(qsub_str)
        qsub.close()

        print("Running " + qsub_path)
        run("{0} {1}".format(cmd, qsub_path), shell=True)
        print("Finished running")
