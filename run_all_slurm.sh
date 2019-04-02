#!/usr/bin/env bash
for EXPERIMENT in ./Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  screen -dmS ${NAME} bash -c "cd ${EXPERIMENT} && ./snakeslurm.sh"
done