import os
from glob import glob
from plot_module import *

for simulated in [True, False]:
    os.chdir("../Data" + ("Simulated" if simulated else "Empirical"))
    exp_dirs = sorted(
        os.listdir("Experiments") if simulated else set(["_".join(i.split("_")[:-1]) for i in os.listdir("Experiments")]))

    for exp in exp_dirs:
        data = dict()
        if simulated:
            for rep, exp_replicate in enumerate(sorted(glob("Experiments/" + exp + "/inference*traces"))):
                for chain, tsv_path in enumerate(sorted(glob(exp_replicate + "/correlation_*.tsv"))):
                    name = "{0}".format(os.path.basename(exp_replicate).split("_")[-2])
                    df = pd.read_csv(tsv_path, sep="\t")
                    if len(df["r2"]) == 0: continue
                    data[name] = df["r2"].values
        else:
            for rep, exp_replicate in enumerate(sorted(glob("Experiments/" + exp + "_*"))):
                for chain, tsv_path in enumerate(sorted(glob(exp_replicate + "/inference*traces/correlation_*.tsv"))):
                    name = "{0} - {1}".format(rep + 1, chain + 1)
                    df = pd.read_csv(tsv_path, sep="\t")
                    if len(df["r2"]) == 0:continue
                    data[name] = df["r2"].values

        if len(data) == 0: continue
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 2.5))
        # plot violin plot
        print(exp)
        for k, v in data.items():
            print("Mean of {0} for {1}".format(v.mean(), k))
        ax.violinplot(data.values(), showmeans=True, showmedians=False)
        ax.set_title('Distribution of $r^2$ between $\\mu$ and $N_{\\mathrm{e}}$ across branches of the tree')
        ax.yaxis.grid(True)
        ax.set_xticks([y + 1 for y in range(len(data))])
        ax.set_xlabel('Model' if simulated else 'Replicate - Chain')
        plt.setp(ax, xticks=[y + 1 for y in range(len(data))],
                 xticklabels=list(data.keys()))
        plt.ylim((0, 1.0))
        ax.set_ylabel('$r^2$ between $\\mu$ and $N_{\\mathrm{e}}$')
        plt.tight_layout()
        plt.savefig("Analysis/identifiability.{0}.pdf".format(exp), format="pdf")
        plt.clf()
