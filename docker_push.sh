#!/usr/bin/env bash
IMAGE_NAME=mutation-selection-drift
TAG=beta
REPO=thibaultlatrille/${IMAGE_NAME}:${TAG}

echo "Building ${REPO}"
docker build -t ${REPO} -f Dockerfile .

read -p "Do you want to push the image? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    docker push ${REPO}
fi