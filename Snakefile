import os
import numpy as np
import pandas as pd
import itertools
import sys

DATASIMULATED = os.path.abspath('./DataSimulated')
np.random.seed(seed=0)
REPLICATES = np.random.choice(9999, 200, replace=False)
os.makedirs(DATASIMULATED + "/matrices", exist_ok=True)

rule all:
    input: expand(DATASIMULATED + '/Experiments/tr20_df1_{rep}', rep=REPLICATES)

rule make_matrix:
    output: DATASIMULATED + '/matrices/tr20_df1_{rep}.tsv'
    input: 'config.yaml'
    run:
        np.random.seed(seed=int(wildcards.rep))
        kappa = 1.0
        dimension = 3
        dof = 1
        sigma_0 = kappa * np.identity(dimension)
        draws = np.random.multivariate_normal(np.zeros(dimension), np.linalg.inv(sigma_0), size=dimension + dof)
        m = np.zeros((dimension, dimension))
        for d in draws:
            d = np.array([d])
            m += np.dot(d.T, d)

        out = dict()
        for i in range(dimension):
            for j in range(i + 1):
                out["c_{1}{0}".format(i, j)] = ["{0:.2g}".format(m[i, j])]

        pd.DataFrame(out).to_csv(str(output), index=False, sep="\t")

rule run_experiment:
    output: directory(DATASIMULATED + '/Experiments/tr20_df1_{rep}')
    input: DATASIMULATED + '/matrices/tr20_df1_{rep}.tsv'
    params: rep=lambda wildcards: str(wildcards.rep),
    shell: 'cp config.yaml DataSimulated && cd DataSimulated && sed -i "s#SEED: 0#SEED: {params.rep}#g" config.yaml && sed -i "s#PRECISION_MATRIX: 0#PRECISION_MATRIX: DataSimulated/matrices/tr20_df1_{params.rep}.tsv#g" config.yaml && python3 simulated_experiment.py -e "tr20_df1_{params.rep}" -c 1 && ln -f -s /home/thibault/MutationSelectionDrift/DataSimulated/run_coevol.sh Experiments/tr20_df1_{params.rep} && cd Experiments/tr20_df1_{params.rep} && sh ./run_coevol.sh'

