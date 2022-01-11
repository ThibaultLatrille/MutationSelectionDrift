# Mutation-Selection-Drift

This repository is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript.
The experiments can either run on a local computer or in a cluster configuration (slurm).

Moreover, this repository gives the tools to produce your own experiment on your dataset, given you have at least 
a DNA alignment file and an associated rooted tree topology (branch lengths are not required).

The experiments are meant to run on Linux/Unix/MacOS operating systems.
However a docker container is provided and a WindowsOS installation (not tested by the author) should in principle work too, at your own peril.

If problems and/or questions are encountered, feel free to [open issues](https://github.com/ThibaultLatrille/MutationSelectionDrift/issues).

## 0. Local copy
Clone the repository and cd to the dir.
```
git clone https://github.com/ThibaultLatrille/MutationSelectionDrift
cd MutationSelectionDrift
```

## 1. Installation
At this step, you may have to chose between either a blue or a red pill :pill:.

With the blue pill, simply pull the prebuild docker container and run the experiments inside the docker container.
Access inside the docker container is also provided by a Jupyter server, hence you will be able to set-up the environment in 2 lines of code,
and start experimenting in your browser.

Or chose the red pill and install the requirements and compile the C++ code for *SimuEvol* and *BayesCode* in your Debian based OS.

Blue pill is prefered is you want minimal conflict with your local system.
Red pill method is prefered if you plan to extensively use the program and tinker with the code.
The two pills are mutually not exclusive, no overdose had ever been observed (though no statistical study had been performed).

### 1.a. Blue pill - Docker
A installation of Docker is required, see https://docs.docker.com/install/.

In Ubuntu (>17.10) with snap pre-installed, one can simply install docker with: 
```
sudo snap install docker
```
Once installed, pull the docker container.
```
docker pull thibaultlatrille/mutation-selection-drift
```
And run the docker container with Jupyter Lab server. 
```
docker run -i -t -p 8888:8888 -v $(pwd)/DataEmpirical:/MutationSelectionDrift/DataEmpirical -v $(pwd)/DataSimulated:/MutationSelectionDrift/DataSimulated thibaultlatrille/mutation-selection-drift /bin/bash -c "jupyter lab --ip='*' --port=8888 --no-browser --allow-root"
```
Once this is done, you have access to a terminal and Jupyter Notebooks from your browser.

To note, the local sub-folders *./DataEmpirical* and *./DataSimulated* are kept synced within the docker container, meaning everything you do in those sub-folders inside the docker container are permanent in your system (and vice-versa).
### 1.b. Red pill - Installation on debian
Install the compiling toolchains:
```
sudo apt install -qq -y make cmake clang openmpi-bin openmpi-common libopenmpi-dev
```
Clone and compile the C++ code for *BayesCode*
```
git clone https://github.com/ThibaultLatrille/bayescode && cd bayescode && git checkout v1.0 && make release && cd ..
```
Clone and compile the C++ code for *SimuEvol*
```
git clone https://github.com/ThibaultLatrille/SimuEvol && cd SimuEvol && git checkout v1.0 && make release && cd ..
```
Install python3 packages
```
sudo apt install -qq -y python3-dev python3-pip screen
pip3 install jupyterlab snakemake numpy matplotlib statsmodels pandas ete3 --user
```
If you wish, you can also run a Jupyter Lab to open the notebooks.
```
jupyter lab
```
## 2. Replicate experiments

To replicate figure 2 of the manuscript on simulated dataset, open the jupyter notebook [ReplicateExperiment.ipynb](https://github.com/ThibaultLatrille/MutationSelectionDrift/blob/master/DataSimulated/ReplicateExperiment.ipynb) in the sub-folder *DataSimulated*.

To replicate figure 3 of the manuscript on mammalian dataset, open the jupyter notebook [ReplicateExperiments.ipynb](https://github.com/ThibaultLatrille/MutationSelectionDrift/blob/master/DataEmpirical/ReplicateExperiments.ipynb) in the sub-folder *DataEmpirical*.

## 3. Run your own experiments 

Instructions can be found in the jupyter notebook [ReplicateExperiments.ipynb](https://github.com/ThibaultLatrille/MutationSelectionDrift/blob/master/DataEmpirical/ReplicateExperiments.ipynb) in the sub-folder *DataEmpirical*.

## 4. Add features or debug in the python scripts
You made modifications to one of the python script, a notebook, this README.md, or you added new features.
You wish this work benefits to all (futur) users of this repository?
Please, feel free to open a [pull-request](https://github.com/ThibaultLatrille/MutationSelectionDrift/pulls)

## 5. Add features or debug in *BayesCode* and *SimuEvol*
You made modifications to the C++ code of either the inference framework *BayesCode* or the simulation framework *SimuEvol*.
You wish this changes benefit to all users of these software?

Please, feel free to open pull-requests in their respective GitHub repository:
* https://github.com/ThibaultLatrille/SimuEvol 
* https://github.com/ThibaultLatrille/bayescode

## Licence

The MIT License (MIT)

Copyright (c) 2019 Thibault Latrille

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
