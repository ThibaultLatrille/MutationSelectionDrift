#!/usr/bin/env bash
git reset --hard
git pull
git clone https://github.com/ThibaultLatrille/SimuEvol
cd SimuEvol
git reset --hard
git pull
cd ..
git clone https://github.com/bayesiancook/bayescode
cd bayescode
git checkout chronogram
git reset --hard
git pull
cd ..
