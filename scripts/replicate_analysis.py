import os
from glob import glob
from ete3 import Tree
from plot_module import *

os.chdir("../DataEmpirical")
exp_dirs = sorted(
    set(["_".join(i.split("_")[:-1]) for i in os.listdir("Experiments") if ("Replicates" in i) and ("eOrthoMam" in i)]))
print(exp_dirs)
os.makedirs("Analysis", exist_ok=True)


def remove_units(str):
    return str.replace("(g)", "").replace("(days)", "").replace("(yrs)", "").replace("(kg)", "").replace("(cm)", "")


dataset = pd.read_csv("Isopods/dataset.tsv", delimiter="\t")


for exp in exp_dirs:
    nhx_dict = dict()
    for exp_replicate in sorted(glob("Experiments/" + exp + "_*")):
        for nhx_path in sorted(glob(exp_replicate + "/*run*.nhx")):
            if nhx_path.count(" ") > 0:
                continue
            nhx_file = os.path.basename(nhx_path)
            if nhx_file not in nhx_dict:
                nhx_dict[nhx_file] = list()
            nhx_dict[nhx_file].append(exp_replicate)

    os.makedirs("Analysis/" + exp, exist_ok=True)
    for tree_name, folder_list in nhx_dict.items():
        print(tree_name, folder_list)
        feature = remove_units(tree_name.split(".")[-2])
        axis_dict, err_dict, taxa_dict = dict(), dict(), dict()

        nb_values = -1
        for index, folder in enumerate(folder_list):
            name = folder.split("_")[-1]
            with open(folder + "/" + tree_name, 'r') as tree_file:
                tree = Tree(remove_units(tree_file.readline()), format=1)

            values = np.array([to_float(getattr(n, feature)) for n in tree.traverse() if feature in n.features])
            taxa_dict[name] = np.array([to_float(getattr(n, feature)) for n in tree if feature in n.features])
            if nb_values == -1:
                nb_values = len(values)
            if len(values) != nb_values:
                print(name + " don't have the same number of values")
                continue
            min_values = np.array(
                [to_float(getattr(n, feature + "_min")) for n in tree.traverse() if feature + "_min" in n.features])
            max_values = np.array(
                [to_float(getattr(n, feature + "_max")) for n in tree.traverse() if feature + "_max" in n.features])
            axis_dict[name] = values
            err_dict[name] = np.vstack((np.abs(values - min_values), np.abs(max_values - values)))

        if "Isopods" in exp:
            df = pd.DataFrame.from_dict(taxa_dict)
            for x in ["eco", "pig", "eye"]:
                df[x] = pd.Series([dataset[dataset['code'] == i].iloc[0][x] for i in tree.get_leaf_names()])
                df.boxplot(column=list(taxa_dict.keys()), by=x, figsize=(len(taxa_dict), 2), layout=(1, len(taxa_dict)), notch=True)
                plt.tight_layout()
                plt.savefig("Analysis/{0}/{1}.{2}.svg".format(exp, tree_name, x))
                plt.clf()
        plot_correlation("Analysis/" + exp + "/" + tree_name.replace(".nhx", ".svg"), axis_dict, err_dict, [],
                         global_xy=True)
