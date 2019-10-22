from Bio import SeqIO
from Bio import Phylo

path = "../DataEmpirical/Vertebrate_25sp"
tree = "rootedtree.nhx"
fasta = "concatenatCDS.fa"

species = [s.name for s in Phylo.read(path + "/" + tree, "newick").get_terminals()]

fasta_file = path + "/" + fasta
fasta_seqs = SeqIO.parse(open(fasta_file, 'r'), 'fasta')

ids_seqs = [(fasta.id, str(fasta.seq)) for fasta in fasta_seqs if (fasta.id in species)]

ali_file = open(fasta_file[:fasta_file.rfind(".")] + ".ali", 'w')
ali_file.write(str(len(ids_seqs)) + " " + str(len(ids_seqs[0][1])) + "\n")
ali_file.write("\n".join([" ".join(id_seq) for id_seq in ids_seqs]))
ali_file.close()

print('Job completed')
