#!/usr/bin/env bash
for EXPERIMENT in ./Experiments/plac_coe*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --touch
  rm -rf all_trace
  rm -rf inference_plot_SimuDiv
  snakemake --printshellcmds -j 4 all_trace
  cd ../..
done