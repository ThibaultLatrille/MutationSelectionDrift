#!/usr/bin/env bash
for EXPERIMENT in ./Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT} && snakemake --printshellcmds
done