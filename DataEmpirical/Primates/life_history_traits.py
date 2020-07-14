#!python3
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt


def transform(x):
    try:
        return np.log(float(x))
    except ValueError:
        return x


traits = pd.read_csv("life_history_traits_exp.tsv", sep="\t")
log_traits = traits.applymap(transform)
log_traits.to_csv("life_history_traits.tsv", index=False, na_rep="NaN", sep="\t")

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
