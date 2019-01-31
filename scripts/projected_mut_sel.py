'''  
T. Latrille (Thibault.latrille@ens-lyon.org), adapted from SJS (stephanie.spielman@gmail.com, https://github.com/clauswilke/Omega_MutSel).
Script to generate Hyphy batchfile, for a variety of model parameterizations, from an alignment in .fasta and a tree.

USAGE: Usage: python prefs_to_freqs.py <fasta> <newick>. The first argument is the alignment file. The second argument is the newick tree file

Dependencies - numpy, jinja2


What does the code do?
First, we find the F61, F1x4, F3x4 frequencies. We use the global alignment frequencies.
Second, we set up the hyphy batch file which makes use of these frequencies.
Third, we generate the MG1 and MG3 matrix files.
'''
from analysis import codon_table, codons, nucindex
import jinja2
import argparse
from ete3 import Tree


def get_nuc_diff(source, target, grab_position=False):
    diff = ''
    position = 5
    for i in range(3):
        if source[i] != target[i]:
            diff += source[i] + target[i]
            position = i
    if grab_position:
        return diff, position
    else:
        return diff


def array_to_hyphy_freq(f):
    ''' Convert array of codon frequencies to a hyphy frequency string to directly write into a batchfile. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{" + str(freq) + "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f


def build_nuc_vars(nuc_freqs, vars_list, constrains_list):
    ''' Compute codon frequencies from GTR nucleotide frequencies. '''
    const_freq = nuc_freqs['T']
    values_set = set(nuc_freqs.values())
    assert len(values_set) == 2 or len(values_set) == 4, "There must be either 2 or 4 frequencies parameters"
    const_freq_sum = "+".join([var for var in values_set if var != const_freq])
    ratio = len(values_set) / len(nuc_freqs.values())

    for var in values_set:
        if var != const_freq:
            vars_list.append("global {0}=0.25;".format(var))

    constrains_list.append("global {0}:={1}-({2});".format(const_freq, ratio, const_freq_sum))
    return 1


def weak_strong(nuc):
    if nuc == "A" or nuc == "T":
        return "W"
    else:
        return "S"


def is_TI(source, target):
    purines = ["A", "G"]
    pyrims = ["C", "T"]
    check1 = source in purines and target in purines
    check2 = source in pyrims and target in pyrims
    if check1 or check2:
        return True
    else:
        return False


def build_rates(param, vars_list, constrains_list):
    vars_set = set()
    assert param in [0, 1, 5]
    gtr_vars = {}
    for n_1 in nucindex.keys():
        for n_2 in nucindex.keys():
            if n_1 != n_2:
                key = n_1 + n_2
                if param == 1 and is_TI(n_1, n_2):
                    gtr_vars[key] = "k"
                if param == 5:
                    value = "exch" + "".join(sorted(n_1 + n_2))
                    if value != 'exchAC':
                        gtr_vars[key] = value
                        vars_set.add(value)
    vars_list.extend(["global {0}=1.0;".format(v) for v in vars_set])
    return gtr_vars


def epsilon_name(codon, omega_param):
    epsilon = "eps" + codon_table[codon]
    if codon_table[codon] == codon_table["ATG"] and omega_param == 95:
        epsilon += "CST"
    return epsilon


def build_matrices(nuc_freqs, exchan_vars, omega_param, vars_list, constrains_list):
    ''' Create MG94-style matrices (use target nucleotide frequencies).  '''
    matrix = '{61, 61, \n'  # MG94
    codon_freqs = [""] * 61
    epsilon_set = set()
    omega_set = set()

    for i, source in enumerate(codons):
        codon_freqs[i] = "*".join([nuc_freqs[source[j]] for j in range(3)])
        if omega_param == 95 or omega_param == 20:
            epsilon = epsilon_name(source, omega_param)
            codon_freqs[i] += "*" + epsilon
            epsilon_set.add(epsilon)

        for j, target in enumerate(codons):

            diff, x = get_nuc_diff(source, target, grab_position=True)
            if len(diff) == 2:
                assert (len(str(x)) == 1), "Problem with determining nucleotide difference between codons"
                freq = str(nuc_freqs[diff[1]])
                # Create string for matrix element
                element = '{' + str(i) + ',' + str(j) + ',mu*t'
                if diff in exchan_vars:
                    element += '*' + exchan_vars[diff]
                if codon_table[source] != codon_table[target]:
                    if omega_param == 1:
                        omega = 'w'
                    elif omega_param == 4:
                        omega = "w_"
                        for nuc in diff:
                            omega += weak_strong(nuc)
                    elif omega_param == 95 or omega_param == 20:
                        if omega_param == 95:
                            omega = 'b_' + "".join(sorted(codon_table[source] + codon_table[target]))
                        epsilon = epsilon_name(target, omega_param)
                        element += '*' + epsilon
                    omega_set.add(omega)
                    element += '*' + omega

                matrix += element + '*' + freq + '} '
    matrix += '}\n'

    if omega_param == 95 or omega_param == 20:
        if omega_param == 95:
            const_epsilon = epsilon_name("ATG", omega_param)
            epsilon_set.remove(const_epsilon)
            assert len(epsilon_set) == 19, "There must be 20 amino-acids"
            const_freq_sum = "+".join([var for var in epsilon_set])
            constrains_list.append("global {0}:=20-({1});".format(const_epsilon, const_freq_sum))

        vars_list.extend(["global {0}=1.0;".format(v) for v in epsilon_set])

    vars_list.extend(["global {0}=1.0;".format(v) for v in omega_set])

    print("{0} omega parameters out of {1:.0f} possible".format(len(omega_set), 20 * 19 / 2))
    constrains_list.append("global z:={0};".format("+".join(codon_freqs)))

    for i, codon in enumerate(codon_freqs):
        codon_freqs[i] = codon + "/z"
    # And save to file.
    return matrix, array_to_hyphy_freq(codon_freqs)


def build_hyphy_batchfile(raw_batch_path, batch_outfile, fasta_infile, tree_infile,
                          rate_param=0, freq_param=3, omega_param=1):
    # Parse input arguments and set up input/outfile files accordingly
    name = batch_outfile.split('/')[-1]
    raw_batchfile = raw_batch_path + "/projected_mut_sel.bf"

    # Calculate frequency parameterizations
    print("Calculating frequency parametrization.")
    nuc_freqs_dict = dict()
    nuc_freqs_dict[1] = {'A': 'pnAT', 'C': 'pnCG', 'G': 'pnCG', 'T': 'pnAT'}
    nuc_freqs_dict[3] = {'A': 'pnA', 'C': 'pnC', 'G': 'pnG', 'T': 'pnT'}

    vars_list = list()
    constrains_list = list()

    vars_list.append("global mu=1;")

    nuc_freqs = nuc_freqs_dict[freq_param]
    build_nuc_vars(nuc_freqs_dict[freq_param], vars_list, constrains_list)

    exchan_vars = build_rates(rate_param, vars_list, constrains_list)

    matrix, codon_freqs = build_matrices(nuc_freqs, exchan_vars, omega_param, vars_list, constrains_list)

    # Create the hyphy batchfile to include the frequencies calculated here. Note that we need to do this since no
    # actual alignment exists which includes all protein positions, so cannot be read in by hyphy file.
    tree_file = Tree(tree_infile)
    tree = tree_file.write(format=0)

    env = jinja2.Environment(loader=jinja2.FileSystemLoader("/"))
    template = env.get_template(raw_batchfile)
    template.stream(matrix=matrix,
                    codon_freqs=codon_freqs,
                    param="r{0}_f{1}_w{2}".format(rate_param, freq_param, omega_param),
                    vars_list=vars_list,
                    constrains_list=constrains_list,
                    tree=tree, name=name,
                    fasta_infile=fasta_infile).dump(batch_outfile)

    print("Building and saving Hyphy Batchfile ({0})".format(batch_outfile))
    return batch_outfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--directory', required=False, type=str,
                        default='/home/thibault/SimuEvol/scripts',
                        dest="d", metavar="<dir_raw_batch>",
                        help="The path to the directory containing the raw batch (projected_mut_sel.bf)")
    parser.add_argument('-b', '--batch', required=False, type=str,
                        default='/home/thibault/SimuEvol/data_hyphy/npcst.bf',
                        dest="b", metavar="<batch_output_path.bf>",
                        help="The path to the output batch.")
    parser.add_argument('-f', '--fasta', required=False, type=str,
                        default="/home/thibault/SimuEvol/data_alignment/npcst.fasta",
                        dest="f", metavar="<.fasta>",
                        help="The path to fasta alignment file")
    parser.add_argument('-t', '--tree', required=False, type=str,
                        default="/home/thibault/SimuEvol/data_trees/np.newick",
                        dest="t", metavar="<.newick>",
                        help="The path to the newick tree file")
    parser.add_argument('-p', '--parameters', required=False, type=str,
                        default="0-3-1",
                        dest="p", metavar="<parameters>",
                        help="The parameters of the GTR-Matrix")
    args = parser.parse_args()
    params = [int(p) for p in args.p.split("-")]
    build_hyphy_batchfile(args.d, args.b, args.f, args.t, params[0], params[1], params[2])
