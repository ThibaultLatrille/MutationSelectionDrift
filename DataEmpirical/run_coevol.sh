#!/usr/bin/env bash
mkdir -p coevol_df0/traces
mkdir -p coevol_df1/traces
mkdir -p coevol_df2/traces
/home/thibault/MutationSelectionDrift/coevol/data/ancov -f -df 0 -t rootedtree.nwk.scaled.nwk -c life_history_traits_exp.traits -x 1 3000 coevol_df0/chain & /home/thibault/MutationSelectionDrift/coevol/data/ancov -f -df 1 -t rootedtree.nwk.scaled.nwk -c life_history_traits_exp.traits -x 1 3000 coevol_df1/chain & /home/thibault/MutationSelectionDrift/coevol/data/ancov -f -df 2 -t rootedtree.nwk.scaled.nwk -c life_history_traits_exp.traits -x 1 3000 coevol_df2/chain &
wait
/home/thibault/MutationSelectionDrift/coevol/data/readancov -x 1500 1 -1 +log coevol_df0/chain & /home/thibault/MutationSelectionDrift/coevol/data/readancov -x 1500 1 -1 +log coevol_df1/chain & /home/thibault/MutationSelectionDrift/coevol/data/readancov -x 1500 1 -1 +log coevol_df2/chain &
wait
python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace coevol_df0/chain.trace --output coevol_df0/traces & python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace coevol_df1/chain.trace --output coevol_df1/traces & python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace coevol_df2/chain.trace --output coevol_df2/traces &
wait
