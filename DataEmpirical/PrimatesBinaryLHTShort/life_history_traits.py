#!python3
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from ete3 import Tree


def t_log(x):
    try:
        return np.log(float(x))
    except ValueError:
        return x


traits = pd.read_csv("lht_exp.tsv", sep="\t")
cols = list(traits)[1:]
filtered_rows = []
for i, row in traits.iterrows():
    if not row[cols].isnull().values.any():
        filtered_rows.append(row)
    else:
        print(row["TaxonName"] + " has not complete data.")
        print(row[cols])
    print()

traits = pd.DataFrame(filtered_rows, columns=list(traits))
traits.to_csv("life_history_traits_exp.tsv", index=False, na_rep="NaN", sep="\t")
print("{0} species with at least one life history trait saved in file 'lht_exp.tsv'.".format(
    len(filtered_rows)))

file = open("life_history_traits_exp.traits", 'w')
file.write("#TRAITS\n{0} {1} ".format(len(traits), len(list(traits)) - 1) +
           " ".join([i for i in list(traits)[1:]]) + "\n" +
           traits.to_csv(index=False, na_rep="-1", sep=" ", header=False))
file.close()
log_traits = traits.applymap(t_log)
log_traits.to_csv("life_history_traits.tsv", index=False, na_rep="NaN", sep="\t")

tree = Tree("full.rootedtree.nwk", format=3)
tree.prune(traits["TaxonName"], preserve_branch_length=True)
tree.write(outfile="rootedtree.nwk", format=3)
pd.DataFrame(traits["TaxonName"]).to_csv("taxlist.lht.csv", index=False, header=None)
print("{0} species with complete data saved in tree 'rootedtree.lht.nwk'.".format(len(traits["TaxonName"])))

axis_dict = {n: log_traits[n].values for n in log_traits if n != "TaxonName"}

my_dpi = 64
f, axs = plt.subplots(len(axis_dict), len(axis_dict), figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)

for row, (row_filename, row_axis) in enumerate(axis_dict.items()):
    for col, (col_filename, col_axis) in enumerate(axis_dict.items()):
        ax = axs[row][col]
        f = np.isfinite(row_axis) & np.isfinite(col_axis)
        ax.scatter(col_axis[f], row_axis[f], label=r"${0}$ points".format(sum(f)))
        idf = np.linspace(min(col_axis[f]), max(col_axis[f]), 100)

        if row_filename != col_filename:
            model = sm.OLS(row_axis[f], sm.add_constant(col_axis[f]))
            results = model.fit()
            b, a = results.params[0:2]
            ax.plot(idf, a * idf + b, '-',
                    label=r"$y={0:.2g}x {3} {1:.2g}$ ($r^2={2:.2g})$".format(a, abs(b), results.rsquared,
                                                                             "+" if float(b) > 0 else "-"))
            ax.legend()
        if row == len(axis_dict) - 1:
            ax.set_xlabel(col_filename)
        if col == 0:
            ax.set_ylabel(row_filename)

plt.savefig("lht.correlation.pdf", format="pdf")
