#!python3
import pandas as pd
import os
import re


def save_df(name, col):
    pd.DataFrame(col).to_csv(name, index=False, header=None)
    print("{0} CDS saved into '{1}'".format(len(col), name))


folder_path = os.path.abspath('..') + "/DataEmpirical/Cetacea"
assert (os.path.isdir(folder_path))
txt = "".join(open("{0}/10_genes_SortaDate_10parts.phy".format(folder_path), "r").readlines())
aligns = txt.split("\n\n")
list_folder_cds = []
for part, align in enumerate(aligns):
    txt = re.sub("\s\s+", " ", align.replace("\n \n", "\n").strip())
    name = "{0}_part_{1}".format(len(aligns), part)
    ali_file = open("{0}/singlegene_alignments/{1}.ali".format(folder_path, name), 'w')
    ali_file.write(txt)
    ali_file.close()
    list_folder_cds.append(name)

print("{0} CDS in the folder 'singlegene_alignments'.".format(len(list_folder_cds)))
save_df("{0}/cds.{1}_part.list".format(folder_path, len(aligns)), list_folder_cds)
