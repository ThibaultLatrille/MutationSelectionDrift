from glob import glob
from csv import reader
import matplotlib.pyplot as plt
import numpy as np

current_dir = "/home/thibault/SimuEvol"
cmd = "sh"
dict_values = dict()

output_path = current_dir + "/polymorphism.tsv"
output = open(output_path, "w")
output.write("Run\tNe\tu\tTheta\tMean(Theta)\tStd(Theta)\tPrecision(%)\tNumberPoints\n")
for folder_path in sorted(glob("{0}/polymorphism_*".format(current_dir))):
    protein = folder_path.split("_")[-3]
    pop_size = int(folder_path.split("_")[-1])

    trace_file_list = sorted(glob("{0}/*.trace".format(folder_path)))
    if len(trace_file_list) > 0:
        for trace_path in trace_file_list:
            name = protein + "_" + str(pop_size) + "_" + trace_path.replace(".trace", "").split("_")[-1]

            with open(trace_path, 'r') as tracefile:
                tracein = list(reader(tracefile, delimiter='\t'))

                for index, col in enumerate(tracein[0]):
                    if col not in dict_values:
                        dict_values[col] = dict()

                    dict_values[col][name] = [float(line[index]) for line in tracein[1:] if float(line[index]) > -6000]

            if len(dict_values["theta"][name]) > 0:
                points = len(dict_values["theta"][name])
                mean_theta = np.mean(dict_values["theta"][name])
                std_theta = 1.96 * np.std(dict_values["theta"][name])
                theta = 4 * pop_size * 1e-6
                output.write(name + "\t" + "\t".join(
                    ["{0:.2e}".format(v) for v in [pop_size, 1e-6, theta, mean_theta, std_theta,
                                                   100 * abs(mean_theta - theta) / theta]]) + "\t" + str(points) + "\n")

output.close()
for col, values in dict_values.items():
    my_dpi = 92
    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    for label, y in values.items():
        plt.plot(range(1, len(y) + 1), y, label=label)
    plt.xlabel('point', fontsize=24)
    plt.ylabel(col, fontsize=24)
    plt.title(col, fontsize=24)
    plt.legend(fontsize=24)
    plt.show()

print("Finished running")
