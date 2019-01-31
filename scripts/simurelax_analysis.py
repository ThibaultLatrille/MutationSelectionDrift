from glob import glob
from csv import reader
from math import isfinite
from numpy import nanmean, nanstd
import matplotlib.pyplot as plt

folder = "data_relax"
folder_path = "/home/thibault/SimuEvol/{0}".format(folder)

for tsv_path in sorted(glob("{0}/*.tsv".format(folder_path))):
    with open(tsv_path, 'r') as tsvfile:
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
        title = "$Ne={0}$, $k={1}$, $\\mu={2}$, $r={3}$, $n={4}$, $m={5}$, $a={6}$, $Q={7}$"
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
                maximum = 1.0
                for label, dict_data in dict_label_data.items():
                    if index:
                        index = len([t for t in dict_data["x"] if t < time])
                    plt.plot(dict_data["x"][index:], dict_data["y"][index:], label=label)
                    maximum = max(max(dict_data["y"][index:]), maximum)
                    print("{0}={1:.3f}±{2:.3f}".format(label, nanmean(dict_data["y"][index:]),
                                                       nanstd(dict_data["y"][index:])))
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                plt.ylim((0, maximum))
                plt.xlabel('$t$', fontsize=16)
                plt.legend()

        plt.suptitle(title)
        save_dir = folder_path + "/" + tsv_path.split("/")[-1].replace("tsv", "png")
        plt.savefig(save_dir, format="png")
        # plt.show()
