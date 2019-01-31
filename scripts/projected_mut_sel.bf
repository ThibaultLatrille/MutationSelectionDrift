/* T. Latrille.
Hyphy inference for an "experimental" dataset using GTR mutation matrix
*/
{% for var in vars_list %}{{ var }} {% endfor %}
{% for var in constrains_list %}{{ var }} {% endfor %}

LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
OPTIMIZATION_PRECSION = 0.00000001;

Matrix_codons={{ matrix }};

/* Read in the data */
DataSet	raw_data = ReadDataFile("{{ fasta_infile }}");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");


/* Set up frequencies. Note that these were all hard-coded in when the file was created via the script SimuEvol/scripts/prefs_to_freqs.py */

Freq_codons={{ codon_freqs }};

/* Optimize likelihoods for each frequency specification */

////////////// {{ param }} //////////////
{% for var in vars_list %}{{ var }} {% endfor %}
{% for var in constrains_list %}{{ var }} {% endfor %}

Model Model_codons = (Matrix_codons, Freq_codons, 0);
UseModel (USE_NO_MODEL);
UseModel(Model_codons);
Tree    Tree_newick = {{ tree }};

branchNames 	= BranchName (Tree_newick, -1);
branchLengths	= BranchLength (Tree_newick, -1);

for (k = 0; k < Columns(branchNames)-1; k=k+1)
{
	ExecuteCommands("Tree_newick." + branchNames[k] + ".t:=" + branchLengths[k] + ";");
}

LikelihoodFunction  LikModel_codons = (filt_data, Tree_newick);
Optimize (paramValues, LikModel_codons);
fprintf (stdout, LikModel_codons);
fprintf ("{{ name }}_hyout.txt", LikModel_codons);
