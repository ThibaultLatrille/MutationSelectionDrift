#!/usr/bin/env bash
for EXPERIMENT in ./DataEmpirical/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --touch
  rm -rf *_traces
  snakemake --printshellcmds --rerun-incomplete -j 8
  cd ../../..
done

for EXPERIMENT in ./DataSimulated/Experiments/Gen5*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --touch
  rm -rf *_traces
  snakemake --printshellcmds --rerun-incomplete -j 8
  cd ../../..
done


for ID in {2540527..2540813..1}; do
  scancel $ID
done