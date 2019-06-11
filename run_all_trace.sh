#!/usr/bin/env bash
for EXPERIMENT in ./DataEmpirical/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --touch inference
  snakemake --printshellcmds -j 6
  cd ../../..
done
for EXPERIMENT in ./DataSimulated/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --touch inference
  snakemake --printshellcmds -j 6
  cd ../../..
done