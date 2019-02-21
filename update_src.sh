#!/usr/bin/env bash
git clone https://github.com/ThibaultLatrille/SimuEvol
cd SimuEvol
git pull
cd ..
git clone https://github.com/bayesiancook/bayescode
cd bayescode
git checkout chronogram
git pull
cd ..
