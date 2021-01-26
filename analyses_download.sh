#!/usr/bin/env bash
# The empirical and simulated analyses can take several weeks to run.
# The figshare repository contains the output of the inference model (BayesCode) and analyses on both simulated (SimuEvol) and empirical dataset.
# https://doi.org/10.6084/m9.figshare.13644110.v1

# Simulated experiments
wget https://ndownloader.figshare.com/files/26193713 -O DataSimulated/Experiments.tar.xz
cd DataSimulated && tar -xf Experiments.tar.xz && cd ..

# Empirical experiments
wget https://ndownloader.figshare.com/files/26193701 -O DataEmpirical/Experiments.tar.xz
cd DataEmpirical && tar -xf Experiments.tar.xz && cd ..

# Empirical replicate analyses
wget https://ndownloader.figshare.com/files/26193695 -O DataEmpirical/Analysis/Analysis.tar.xz
cd DataEmpirical/Analysis && tar -xf Analysis.tar.xz && cd ../..


