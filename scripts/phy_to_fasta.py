import os

current_dir = "/home/thibault/SimuEvol/"
folder = "data_primates"
folder_path = "{0}/{1}".format(current_dir, folder)
phy_files = [f for f in os.listdir(folder_path) if ".phy" == f.strip()[-4:]]
print("Found {0} fasta files".format(len(phy_files)))
hit = 0

for phy_file in sorted(phy_files):
    hit += 1
    seq_len = -1
    phy_path = "{0}/{1}".format(folder_path, phy_file)
    fasta_file = open(phy_path.replace(".phy", ".fasta"), 'w')
    read_phy = open(phy_path, 'r')
    read_phy.readline()
    species = 0
    for line_read_phy in read_phy:
        assert line_read_phy.count("\t") == 1
        split_line_read_phy = line_read_phy.split("\t")
        name, seq = split_line_read_phy
        seq = seq.replace("\n", "")
        if seq_len == -1:
            seq_len = len(seq)
            assert seq_len % 3 == 0
        else:
            assert seq_len == len(seq)
        counter = 0
        unknown_chars = seq
        for letter in "ATGC-":
            counter += seq.count(letter)
            unknown_chars = unknown_chars.replace(letter, "")
        if seq_len != counter:
            print("UNKNOWN characters : \n" + unknown_chars)
        fasta_file.write(">" + name + "\n" + seq + "\n")
        species += 1
    fasta_file.close()
    print("{0} species for phy file {1}, with {2} nucleotides".format(species, phy_file, seq_len))

print("{0} phy files".format(hit))
print('Job completed')
