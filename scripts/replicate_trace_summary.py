import os
from glob import glob
import pandas as pd
import numpy as np
from subprocess import run

folder = "DataSimulated"
os.chdir("../" + folder)
os.makedirs("Traces/", exist_ok=True)
param = "aaent"
models = ["node", "aa"]

table = open("Traces/{0}.tex".format(param), 'w')
table.writelines("\\documentclass[USLetter,5pt]{article}\n"
                 "\\usepackage{adjustbox}\n")
table.writelines("\\newcommand{\\specialcell}[2][c]{%\n\\begin{tabular}[#1]{@{}c@{}}#2\\end{tabular}}\n")
table.writelines("\\begin{document}\n")
heading = [param] + models
table.writelines("\\begin{table}[ht]\n\\centering\n\\begin{adjustbox}{width = 1\\textwidth}\n")
table.writelines("\\begin{tabular}{|" + "c|" * len(heading) + "}\n")
table.writelines("\\hline\n")
table.writelines(" & ".join(heading) + "\\\\\n")
table.writelines("\\hline\n")

for exp in sorted([d for d in os.listdir("Experiments")]):
    name_traces = dict()
    for filepath in sorted(glob("Experiments/" + exp + "/*.trace.tsv")):
        models_counts = [m for m in models if os.path.basename(filepath).count(m + "_") > 0]
        if len(models_counts) == 0:
            continue
        else:
            assert (len(models_counts) == 1)
            model = models_counts[0]

        row_index = os.path.basename(filepath).replace(model + "_", "").replace(".trace.tsv", "").replace("_", " ")
        if row_index not in name_traces:
            name_traces[row_index] = dict()
        name_traces[row_index][model] = pd.read_csv(filepath, sep='\t')[param].values

    for name, traces in name_traces.items():
        elts = [exp.replace("_", " ") + " " + name] + [
            "${0:.2f} \\pm {1:.2f}$".format(np.nanmean(traces[m]), 1.96 * np.nanstd(traces[m])) if (m in traces) else ""
            for m in models]
        table.writelines(" & ".join(elts) + "\\\\\n")
        table.writelines("\\hline\n")

table.writelines("\\end{tabular}\n")
table.writelines("\\end{adjustbox}\n" +
                 "\\caption{}\n" + "\\end{table}\n")
table.writelines("\\end{document}\n")
table.close()
run("pdflatex -output-directory {0} {1}.tex".format("Traces", param), shell=True)
