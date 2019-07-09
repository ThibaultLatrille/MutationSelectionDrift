#!python3
import pandas as pd
from ete3 import Tree
import numpy as np
import os

excel_file = pd.ExcelFile("Species.xls")
sheet_names = excel_file.sheet_names
print("{0} species in Species.xls.".format(len(sheet_names)))
cols = ["Longevity (days)", "Length (cm)", "Weight (kg)"]

for folder in ["47SP", "Vertebrates"]:
    print(folder)
    folder_path = os.path.abspath('..') + "/" + folder
    assert (os.path.isdir(folder_path))

    tree = Tree(folder_path + "/rootedtree.nhx", format=3)
    leave_names = tree.get_leaf_names()

    idx = False
    filtered_rows = []
    filtered_tax_list = []
    for taxon in sheet_names:
        if taxon in leave_names:
            df = excel_file.parse(taxon)
            row = df.loc[0, cols]
            if not row.isnull().values.any():
                filtered_tax_list.append(taxon)
            if not row.isnull().values.all():
                filtered_rows.append([taxon] + [np.log(i) for i in row])

    pd.DataFrame(filtered_rows, columns=(["TaxonName"] + [i.replace(" ", "_") for i in cols])).to_csv(
        folder_path + "/life_history_traits.tsv", na_rep="NaN", sep='\t', index=False)
    print(
        "{0} species (out of {1}) with at least one life history trait saved in file 'life_history_traits.tsv'.".format(
            len(filtered_rows), len(leave_names)))

    tree.prune(filtered_tax_list, preserve_branch_length=True)
    tree.write(outfile=folder_path + "/rootedtree.lht.nhx", format=3)
    pd.DataFrame(filtered_tax_list).to_csv(folder_path + "/taxlist.lht.csv", index=False, header=None)
    print(
        "{0} species (out of {1}) with complete data saved in tree 'rootedtree.lht.nhx'.".format(len(filtered_tax_list),
                                                                                                 len(leave_names)))
