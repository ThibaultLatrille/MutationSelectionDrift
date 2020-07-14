#!/usr/bin/env bash
mkdir -p bayescode_diag_df0/traces
mkdir -p bayescode_diag_df1/traces
mkdir -p bayescode_diag_df2/traces
/home/thibault/MutationSelectionDrift/bayescode/_build/nodetraitsmutsel  --df 0 -t rootedtree.nwk.annotated -u 4000 --traitsfile life_history_traits.tsv --fossils rootedtree.tsv bayescode_diag_df0/chain & /home/thibault/MutationSelectionDrift/bayescode/_build/nodetraitsmutsel  --df 1 -t rootedtree.nwk.annotated -u 4000 --traitsfile life_history_traits.tsv --fossils rootedtree.tsv bayescode_diag_df1/chain & /home/thibault/MutationSelectionDrift/bayescode/_build/nodetraitsmutsel  --df 2 -t rootedtree.nwk.annotated -u 4000 --traitsfile life_history_traits.tsv --fossils rootedtree.tsv bayescode_diag_df2/chain &
wait
python3 /home/thibault/MutationSelectionDrift/scripts/plot_correlation_matrix.py -b 2000 --trace bayescode_diag_df0/chain --output ./bayescode_diag_df0/chain.cov & python3 /home/thibault/MutationSelectionDrift/scripts/plot_correlation_matrix.py -b 2000 --trace bayescode_diag_df1/chain --output ./bayescode_diag_df1/chain.cov & python3 /home/thibault/MutationSelectionDrift/scripts/plot_correlation_matrix.py -b 2000 --trace bayescode_diag_df2/chain --output ./bayescode_diag_df2/chain.cov &
wait
/home/thibault/MutationSelectionDrift/bayescode/_build/readnodetraitsmutsel --burnin 2000 --trace bayescode_diag_df0/chain & /home/thibault/MutationSelectionDrift/bayescode/_build/readnodetraitsmutsel --burnin 2000 --trace bayescode_diag_df1/chain & /home/thibault/MutationSelectionDrift/bayescode/_build/readnodetraitsmutsel --burnin 2000 --trace bayescode_diag_df2/chain &
wait
python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace bayescode_diag_df0/chain.trace --output bayescode_diag_df0/traces & python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace bayescode_diag_df0/chain.trace --output bayescode_diag_df1/traces & python3 /home/thibault/MutationSelectionDrift/scripts/plot_traces.py --trace bayescode_diag_df0/chain.trace --output bayescode_diag_df2/traces &
wait