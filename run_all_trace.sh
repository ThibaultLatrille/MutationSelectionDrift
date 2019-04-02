#!/usr/bin/env bash
for EXPERIMENT in ./Experiments/Chronogram*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  # cp Snakefile ${EXPERIMENT}
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --touch
  rm -rf all_trace
  rm -rf inference_plot_SimuDiv
  snakemake --printshellcmds -j 4 all_trace
  cd ../..
done