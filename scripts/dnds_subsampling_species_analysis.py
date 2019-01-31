import os
import pandas as pd
import matplotlib.pyplot as plt

site_wise = False

folders = ["siteomega", "mutsel"]
data_path = "/home/thibault/SimuEvol"
results_path = "{0}/data_pb".format(data_path)

output_path = results_path + 'summary.csv'
output_file = open(output_path, 'w')
header = ["SelectionCoefficient", "ProbaPermutation", "NodeIndex",
          "NbrLeaves", "ChainNbr", "Option", "Omega"]
if site_wise:
    header.extend(["Site"])
output_file.write("\t".join(header) + "\n")

for file in os.listdir(results_path):
    file_path = "{0}/{1}".format(results_path, file)

    run_id = file.replace('.predsiteomega', '').replace('.trace', '').split('_')
    assert len(run_id) == 7
    protein, s, p, index, n, opt_name, chain = run_id

    if site_wise:
        if file.endswith(".predsiteomega"):
            table = pd.read_table(file_path, sep="\t", index_col=0)
            for site in table.index.values:
                output_file.write(
                    "\t".join(map(str, [s, p, index, n, chain, opt_name, table.loc[site, "omega"], site])) + "\n")
        if file.endswith(".trace") and ("siteomega" in file):
            table = pd.read_table(file_path, sep="\t")
            for site in [col for col in table.columns.values if col.isdigit()]:
                output_file.write(
                    "\t".join(map(str, [s, p, index, n, chain, opt_name, table[site].mean(), site])) + "\n")
    else:
        if file.endswith(".predsiteomega") or (
                file.endswith(".trace") and (("mutselfreeomega" in file) or ("globalomega" in file))):
            table = pd.read_table(file_path, sep="\t")
            if file.endswith(".predsiteomega") and ("mutselfreeomega" in file):
                opt_name = "Pred" + opt_name
            output_file.write("\t".join(map(str, [s, p, index, n, chain, opt_name, table["omega"].mean()])) + "\n")
        if file.endswith(".trace") and ("siteomega" in file):
            table = pd.read_table(file_path, sep="\t")
            output_file.write("\t".join(map(str, [s, p, index, n, chain, opt_name, table["meanomega"].mean()])) + "\n")

output_file.close()

df = pd.read_table(output_path, sep="\t")
for s in set(df.SelectionCoefficient):
    for p in set(df.ProbaPermutation):
        filtered_df = df[(df.SelectionCoefficient == s) & (df.ProbaPermutation == p)]
        if len(filtered_df) > 0:
            for opt_name in sorted(set(df.Option), reverse=True):
                opt_df = filtered_df[filtered_df.Option == opt_name]
                plt.scatter(opt_df.NbrLeaves, opt_df.Omega, label=opt_name)

            plt.legend()
            plt.xlabel("NbrLeaves")
            plt.ylabel("Omega")
            plt.ylim((0, max(max(filtered_df.Omega), 1)))
            plt.title("S={0}, P={1}".format(s, p))
            plt.show()
