#!python3
import pandas as pd
import os

orthomam = pd.read_csv("./orthomam.tsv", sep="\t")
print("{0} genes with mutsel inference.".format(len(orthomam)))

genes = orthomam.query("pp_m2a < 0.1 and pp_mutsel < 0.1 and pp_mutselfix < 0.1")["gene"]
print("{0} genes with mutsel inference not under positive selection".format(len(genes)))

folder_list = [i.replace(".ali", "") for i in os.listdir('./singlegene_alignments')]
print("{0} genes in the folder 'singlegene_alignments'.".format(len(folder_list)))

ortho_list = genes[genes.isin(folder_list)]
print("{0} genes in the folder 'singlegene_alignments' and not under positive selection.".format(len(ortho_list)))

filename = "./filtered.list"
ortho_list.to_csv("./" + filename, index=False)
print("Filtered list saved into '{0}'".format(filename))
