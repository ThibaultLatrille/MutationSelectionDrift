#!python3
import pandas as pd
import numpy as np


def transform(x):
    try:
        return np.log(float(x))
    except ValueError:
        return x


traits = pd.read_csv("life_history_traits_exp.tsv", sep="\t")
traits.applymap(transform).to_csv("life_history_traits.tsv", index=False, na_rep="NaN", sep="\t")
