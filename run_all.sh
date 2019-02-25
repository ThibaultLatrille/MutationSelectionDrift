#!/usr/bin/env bash
for EXPERIMENT in ./Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  screen -dmS ${NAME} "cd ${EXPERIMENT} && ./snakeslurm.sh"
done