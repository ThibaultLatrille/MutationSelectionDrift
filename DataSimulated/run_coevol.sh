#!/usr/bin/env bash
python3 /home/thibault/MutationSelectionDrift/scripts/exponentiate_lht.py --input SimuDiv_exp.traits.tsv.drawn
python3 /home/thibault/MutationSelectionDrift/scripts/calibs_tree.py --nwk SimuDiv_exp.tree
mkdir -p coevol_Df2/traces
mkdir -p coevol_Df1/traces
mkdir -p coevol/traces
/home/thibault/MutationSelectionDrift/coevol/data/ancov -f -df 2 -t SimuDiv_exp.tree.scaled.nwk -c SimuDiv_exp.exponentiate.traits.tsv.drawn -x 1 1000 coevol_Df2/chain & /home/thibault/MutationSelectionDrift/coevol/data/ancov -f -df 1 -t SimuDiv_exp.tree.scaled.nwk -c SimuDiv_exp.exponentiate.traits.tsv.drawn -x 1 1000 coevol_Df1/chain & /home/thibault/MutationSelectionDrift/coevol/data/ancov -f -t SimuDiv_exp.tree.scaled.nwk -c SimuDiv_exp.exponentiate.traits.tsv.drawn -x 1 1000 coevol/chain &
wait
/home/thibault/MutationSelectionDrift/coevol/data/readancov -x 500 1 -1 +log coevol_Df2/chain & /home/thibault/MutationSelectionDrift/coevol/data/readancov -x 500 1 -1 +log coevol_Df1/chain & /home/thibault/MutationSelectionDrift/coevol/data/readancov -x 500 1 -1 +log coevol/chain &
wait
