# Dockerfile
FROM debian

# install necessary tools
RUN apt-get update
RUN apt-get install -qq -y git make cmake clang openmpi-bin openmpi-common

# pull git repository
RUN git clone https://github.com/ThibaultLatrille/MutationSelectionDrift
WORKDIR MutationSelectionDrift

# install BayesCode
RUN git clone https://github.com/bayesiancook/bayescode && cd bayescode && git checkout chronogram && make release

# install SimuEvol
RUN git clone https://github.com/ThibaultLatrille/SimuEvol && cd SimuEvol && git checkout 3291eb3e6e89b6a72afe5b87bc0b04cbd72c9f0c && make release

# Install python packages
RUN apt-get install -qq -y python3-dev python3-pip screen
RUN pip3 install jupyterlab snakemake numpy matplotlib statsmodels seaborn pandas ete3