#!/usr/bin/env bash
python3 /home/thibault/MutationSelectionDrift/scripts/exponentiate_lht.py --input SimuDiv_exp.traits.tsv
python3 /home/thibault/MutationSelectionDrift/scripts/calibs_tree.py --nwk SimuDiv_exp.tree
mkdir -p coevol_Df1_PriorSigma1/traces
/home/thibault/MutationSelectionDrift/coevol/data/ancov -f -df 1 -priorsigma 1 -t SimuDiv_exp.tree.scaled.nwk -c SimuDiv_exp.exponentiate.traits.tsv -x 1 3000 coevol_Df1_PriorSigma1/chain
/home/thibault/MutationSelectionDrift/coevol/data/readancov -x 1500 1 -1 +log coevol_Df1_PriorSigma1/chain
python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace coevol_Df1_PriorSigma1/chain.trace --output coevol_Df1_PriorSigma1/traces
mkdir -p coevol_PriorSigma1/traces
/home/thibault/MutationSelectionDrift/coevol/data/ancov -f -priorsigma 1 -t SimuDiv_exp.tree.scaled.nwk -c SimuDiv_exp.exponentiate.traits.tsv -x 1 3000 coevol_PriorSigma1/chain
/home/thibault/MutationSelectionDrift/coevol/data/readancov -x 1500 1 -1 +log coevol_PriorSigma1/chain
python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace coevol_PriorSigma1/chain.trace --output coevol_PriorSigma1/traces
mkdir -p coevol_Df1/traces
/home/thibault/MutationSelectionDrift/coevol/data/ancov -f -df 1 -t SimuDiv_exp.tree.scaled.nwk -c SimuDiv_exp.exponentiate.traits.tsv -x 1 3000 coevol_Df1/chain
/home/thibault/MutationSelectionDrift/coevol/data/readancov -x 1500 1 -1 +log coevol_Df1/chain
python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace coevol_Df1/chain.trace --output coevol_Df1/traces
mkdir -p coevol/traces
/home/thibault/MutationSelectionDrift/coevol/data/ancov -f -t SimuDiv_exp.tree.scaled.nwk -c SimuDiv_exp.exponentiate.traits.tsv -x 1 3000 coevol/chain
/home/thibault/MutationSelectionDrift/coevol/data/readancov -x 1500 1 -1 +log coevol/chain
python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace coevol/chain.trace --output coevol/traces