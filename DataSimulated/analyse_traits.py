#!python3
import pandas as pd
import matplotlib
from glob import glob
import numpy as np
import os
matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

my_dpi = 128
plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)

folder = "tr20_df1_"
col = "PriorSigma_0"
p_vals = []
for exp_replicate in sorted(glob("Experiments/{0}*".format(folder))):
    for trace in sorted(glob(exp_replicate + "/*_run.trace")):
        vals = pd.read_csv(trace, sep='\t', usecols=[col])[col].values[50:]
        if len(vals) < 50: continue
        p_vals.append(np.sum(vals >= 1.0) / len(vals))

p_vals = sorted(p_vals)
x = np.linspace(0, 1, len(p_vals))
plt.scatter(x, p_vals)
plt.plot(x, x, linewidth=2, color="black")
plt.xlabel('Rank')
plt.ylabel('$p$')
plt.tight_layout()
plt.savefig('analyse/{0}bayescode_analyse.svg'.format(folder), format='svg')
plt.clf()

col = "diag0"
for model in ["coevol", "coevol_Df1", "coevol_Df2"]:
    p_vals = []
    for exp_replicate in sorted(glob("Experiments/{0}*".format(folder))):
        file = exp_replicate + "/" + model + "/chain.trace"
        if not os.path.isfile(file): continue
        vals = pd.read_csv(file, sep='\t', usecols=[col])[col].values[500:]
        p_vals.append(np.sum(vals >= 1.0) / len(vals))

    p_vals = sorted(p_vals)
    x = np.linspace(0, 1, len(p_vals))
    plt.scatter(x, p_vals)
    plt.plot(x, x, linewidth=2, color="black")
    plt.xlabel('Rank')
    plt.ylabel('$p$')
    plt.tight_layout()
    plt.savefig('analyse/{0}{1}_analyse.svg'.format(folder, model), format='svg')
    plt.clf()
plt.close('all')
