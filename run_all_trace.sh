#!/usr/bin/env bash
for EXPERIMENT in ./DataEmpirical/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  rm -rf inference_*_trees
  rm -rf ./CorrelationMatrices
  rm -rf ./all_correlation_matrix
  rm -rf all_profiles
  snakemake --unlock
  snakemake --touch inference
  snakemake --printshellcmds --rerun-incomplete -j 4
  cd ../../..
done

for EXPERIMENT in ./DataSimulated/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  rm -rf all_profiles
  rm -rf all_traces
  rm -rf all_tress
  rm -rf simulation_*
  snakemake --unlock
  snakemake --touch inference
  snakemake --printshellcmds --rerun-incomplete -j 4
  cd ../../..
done
