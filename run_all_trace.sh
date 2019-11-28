#!/usr/bin/env bash
for EXPERIMENT in ./DataEmpirical/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake build --rerun-incomplete
  snakemake --touch inference
  snakemake --printshellcmds -j 8 --rerun-incomplete
  cd ../../..
done
for EXPERIMENT in ./DataSimulated/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake build
  snakemake --touch inference
  snakemake --printshellcmds --rerun-incomplete -j 4
  cd ../../..
done