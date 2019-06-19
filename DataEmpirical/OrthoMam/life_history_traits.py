#!python3
import pandas as pd
from ete3 import Tree
import numpy as np

traits = pd.read_csv("anage_data.tsv", sep="\t")
print("{0} species in AnAge Database.".format(len(traits)))

tree = Tree("placnr.rootedtree.nhx", format=3)
leave_names = tree.get_leaf_names()

idx = False
cols = ["Maximum longevity (yrs)", "Adult weight (g)", "Female maturity (days)"]
tax_list = pd.read_csv("placnr.taxlist", sep="\t", header=None)[0].values
filtered_rows = []
filtered_tax_list = []
for v in tax_list:
    assert (v in leave_names)
    genus = v.split("_")[0]
    species = v.split("_")[1]
    row = traits.loc[(traits["Genus"] == genus), cols].mean()
    if not row.isnull().values.any():
        print(v + " has complete data (added to file 'placnr.lht.taxlist').")
        filtered_tax_list.append(v)
    if not row.isnull().values.all():
        print(v + " has been added to file 'life_history_traits.tsv'.")
        filtered_rows.append([v] + list(np.log(row)))
    else:
        print(v + " has no data available.")
        print(row)
    print()

pd.DataFrame(filtered_rows, columns=(["TaxonName"] + [i.replace(" ", "_") for i in cols])).to_csv(
    "life_history_traits.tsv", na_rep="NaN", sep='\t', index=False)
print("{0} species with at least one life history trait saved in file 'life_history_traits.tsv'.".format(
    len(filtered_rows)))

tree.prune(filtered_tax_list, preserve_branch_length=True)
tree.write(outfile="placnr.lht.rootedtree.nhx", format=3)
pd.DataFrame(filtered_tax_list).to_csv("placnr.lht.taxlist", index=False, header=None)
print("{0} species with complete data saved in tree 'placnr.lht.rootedtree.nhx'.".format(len(filtered_tax_list)))

orthomam_v9_tree = Tree("orthomam_v9.nhx")
orthomam_v9_names = set(orthomam_v9_tree.get_leaf_names())
short_names = [i for i in tree.get_leaf_names() if i.split("_")[0] in orthomam_v9_names]
tree.prune(short_names, preserve_branch_length=True)
tree.write(outfile="placnr.short.rootedtree.nhx", format=3)
pd.DataFrame(short_names).to_csv("placnr.short.taxlist", index=False, header=None)
print("{0} species with complete data saved in tree 'placnr.short.rootedtree.nhx'.".format(len(short_names)))
