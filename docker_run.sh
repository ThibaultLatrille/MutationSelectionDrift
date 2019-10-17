#!/usr/bin/env bash
IMAGE_NAME=mutation-selection-drift
TAG=beta
REPO=thibaultlatrille/${IMAGE_NAME}:${TAG}

echo "Building ${REPO}"
docker build -t ${REPO} -f Dockerfile .

echo "Running ${REPO}"
docker run -i -t -p 8888:8888 -v $(pwd)/DataEmpirical:/MutationSelectionDrift/DataEmpirical -v $(pwd)/DataSimulated:/MutationSelectionDrift/DataSimulated ${REPO} /bin/bash -c "jupyter lab --ip='*' --port=8888 --no-browser --allow-root"
