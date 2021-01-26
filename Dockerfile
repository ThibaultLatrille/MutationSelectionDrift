# Dockerfile
FROM debian

# install necessary tools
RUN apt-get update
RUN apt-get install -qq -y git make cmake clang openmpi-bin openmpi-common

# pull git repository
RUN git clone https://github.com/ThibaultLatrille/MutationSelectionDrift
WORKDIR MutationSelectionDrift

# install BayesCode
RUN git clone https://github.com/ThibaultLatrille/bayescode && cd bayescode && make release

# install SimuEvol
RUN git clone https://github.com/ThibaultLatrille/SimuEvol && cd SimuEvol && make release

# Install python packages
RUN apt-get install -qq -y python3-dev python3-pip screen
RUN pip3 install jupyterlab snakemake numpy matplotlib statsmodels seaborn pandas ete3
