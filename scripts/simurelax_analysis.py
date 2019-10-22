from glob import glob
from csv import reader
from math import isfinite
from numpy import nanmean, nanstd
import matplotlib.pyplot as plt

folder_path = "/home/thibault/PolyMutSel/DataSimulated/SimuRelax"

for tsv_path in sorted(glob("{0}/*.tsv".format(folder_path))):
    with open(tsv_path, 'r') as tsvfile:
        print(tsv_path)
        tsvin = list(reader(tsvfile, delimiter='\t'))

        params = tsv_path.split("/")[-1].replace('.tsv', '').split("_")
        population_size = int(params[0])
        k = int(params[1])
        mu = float(params[2])
        r = float(params[3])
        n = int(params[4])
        m = int(params[5])
        a = float(params[6])
        q = float(params[7])
        title = "$Ne={0}$, $k={1}$, $\\mu={2}$, $r={3}$, $c={4}$, $p={5}$, $a={6}$, $Q={7}$"
        title = title.format(population_size, k, mu, r, n, m, a, q)

        time = int([row[1] for row in tsvin if row[0] == "t"][0])

        data = [row for row in tsvin if len(row) == 4 and row[1] != "0"]
        dict_num_label_data = {}
        for row in data:
            if not row[1] in dict_num_label_data:
                dict_num_label_data[row[1]] = {}
            if not row[2] in dict_num_label_data[row[1]]:
                dict_num_label_data[row[1]][row[2]] = {"x": [], "y": []}

            if isfinite(float(row[3])):
                dict_num_label_data[row[1]][row[2]]["x"].append(float(row[0]))
                dict_num_label_data[row[1]][row[2]]["y"].append(float(row[3]))

        my_dpi = 92
        rc_size = len(dict_num_label_data)
        plt.figure(figsize=(1920 / my_dpi, 360 * rc_size / my_dpi), dpi=my_dpi)
        for i, (num, dict_label_data) in enumerate(dict_num_label_data.items()):
            for index in [False, True]:
                plt.subplot(rc_size, 2, i * 2 + index + 1)
                maximum = 0.0
                minimum = 0.0
                for label, dict_data in dict_label_data.items():
                    if index:
                        index = len([t for t in dict_data["x"] if t < time])
                    plt.plot(dict_data["x"][index:], dict_data["y"][index:], label=label)
                    maximum = max(max(dict_data["y"][index:]), maximum)
                    minimum = min(min(dict_data["y"][index:]), minimum)
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                plt.ylim((minimum, maximum))
                plt.xlabel('$t$', fontsize=16)
                plt.legend()
            if "ΔF" in label:
                print("ΔF Selection {0}={1:.3f}±{2:.3f}".format(label, nanmean(dict_data["y"][10:index]),
                                                                nanstd(dict_data["y"][10:index])))
                print("ΔF Relax {0}={1:.3f}±{2:.3f}".format(label, nanmean(dict_data["y"][index:]),
                                                            nanstd(dict_data["y"][index:])))
        plt.suptitle(title)
        save_dir = folder_path + "/" + tsv_path.split("/")[-1].replace("tsv", "svg")
        plt.savefig(save_dir, format="svg")
        # plt.show()
        plt.clf()
        plt.close("all")
