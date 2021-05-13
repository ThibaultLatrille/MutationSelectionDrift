import os
from glob import glob
from ete3 import Tree
from plot_module import *

os.chdir("../DataEmpirical")
exp_dirs = sorted(
    set(["_".join(i.split("_")[:-1]) for i in os.listdir("Experiments") if ("Replicates" in i) and ("Isopods" in i)]))

os.makedirs("Analysis/DataFrame", exist_ok=True)
file_format = "pdf"


def is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


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
        for index, folder in enumerate(sorted(folder_list)):
            name = folder.split("_")[-1]
            with open(folder + "/" + tree_name, 'r') as tree_file:
                tree = Tree(remove_units(tree_file.readline()), format=1)

            values = np.array([to_float(getattr(n, feature)) for n in tree.traverse() if feature in n.features])
            taxa_dict["Rep. {0}".format(index + 1)] = np.array(
                [to_float(getattr(n, feature)) for n in tree if feature in n.features])
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

        plot_correlation("Analysis/" + exp + "/" + tree_name.replace(".nhx", '.' + file_format), axis_dict, err_dict,
                         global_min_max=True)
        if ("LogPopulationSize" not in tree_name) and ("LogMutationRate" not in tree_name): continue
        df = pd.DataFrame.from_dict(taxa_dict)
        df = df.apply(np.exp)
        ratio = ["\\textbf{ " + "{0:.3g}".format(max(df[col]) / min(df[col])) + "}" for col in df]
        if "Isopods" in exp:
            label = label_transform(tree_name.split(".")[-2])
            dict_trait = {"eco": "Habitat", "pig": "Pigmentation", "eye": "Ocular structure"}
            dict_trait_def = {"eco": {"epi": "Surface", "hypo": "Underground"},
                              "pig": {"D": "Depigmented", "DP": "Part. dep.", "P": "Pigmented"},
                              "eye": {"A": "Anophthalmia", "M": "Microphthalmia", "O": "Ocular"}}

            for trait, trait_name in dict_trait.items():
                df[trait_name] = pd.Series([dict_trait_def[trait][dataset[dataset['code'] == i].iloc[0][trait]] for i in
                                            tree.get_leaf_names()])
                axes = df.boxplot(column=list(taxa_dict.keys()), by=trait_name, figsize=(len(taxa_dict), 2),
                                  layout=(1, len(taxa_dict)), rot=90, fontsize=6)
                axes[0].set_ylabel(label, fontsize=6)
                for ax in axes:
                    ax.set_xlabel("")
                    ax.set_title(ax.get_title(), {'fontsize': 8})
                plt.suptitle("")
                plt.tight_layout()
                plt.savefig("Analysis/{0}/{1}.{2}.{3}".format(exp, tree_name, trait, file_format), format=file_format)
                plt.clf()
            df["Code"] = pd.Series(tree.get_leaf_names())
            taxa = [dataset[dataset['code'] == i].iloc[0]["morphosp"] for i in tree.get_leaf_names()]
            taxa = ["\\textit{" + " ".join(str(i).replace('nan', '-').split(" ")[:2]) + "}" for i in taxa]
            ratio += ['-'] * 4
            frames = []
            for k in taxa_dict.keys():
                d = {label: df[k], 'Code': df["Code"], "Rep": [k] * len(df[k])}
                for trait_name in dict_trait.values():
                    d[trait_name] = df[trait_name]
                frames.append(pd.DataFrame(d))
            merged = pd.concat(frames)
            for trait, trait_name in dict_trait.items():
                bp = merged.boxplot(column=[label], by=trait_name, fontsize=18, notch=True,
                                    patch_artist=True, return_type="dict")
                [[item.set_linewidth(3) for item in bp[key]['boxes']] for key in bp.keys()]
                [[item.set_linewidth(3) for item in bp[key]['fliers']] for key in bp.keys()]
                [[item.set_linewidth(3) for item in bp[key]['medians']] for key in bp.keys()]
                [[item.set_linewidth(3) for item in bp[key]['whiskers']] for key in bp.keys()]
                [[item.set_linewidth(3) for item in bp[key]['caps']] for key in bp.keys()]

                [[item.set_color(BLUE) for item in bp[key]['boxes']] for key in bp.keys()]
                [[item.set_color("black") for item in bp[key]['fliers']] for key in bp.keys()]
                [[item.set_color(GREEN) for item in bp[key]['medians']] for key in bp.keys()]
                [[item.set_color(BLUE) for item in bp[key]['whiskers']] for key in bp.keys()]
                [[item.set_color("black") for item in bp[key]['caps']] for key in bp.keys()]

                ax = plt.gca()
                ax.set_ylabel(label, fontsize=22)
                ax.set_xlabel("")
                ax.set_title(trait_name, fontsize=26)
                plt.suptitle("")
                plt.tight_layout()
                plt.savefig("Analysis/{0}/{1}.{2}.merged.{3}".format(exp, tree_name, trait, file_format),
                            format=file_format)
                plt.clf()
            merged.to_csv("Analysis/DataFrame/" + exp + "_merged_" + tree_name.replace(".nhx", '.tsv'), index=False,
                          sep="\t")
            merged.to_latex("Analysis/DataFrame/" + exp + "_merged_" + tree_name.replace(".nhx", '.tex'), index=False,
                            float_format="%.3g")
        else:
            taxa = [(("\\textit{" + i.replace('_', ' ') + "}") if ("_" in i) else i) for i in tree.get_leaf_names()]
        df["Taxon"] = pd.Series(taxa)
        df.loc[len(df.index)] = ratio + ["\\textbf{Maximum range}"]
        df.to_csv("Analysis/DataFrame/" + exp + "_" + tree_name.replace(".nhx", '.tsv'), index=False, sep="\t")
        df.to_latex("Analysis/DataFrame/" + exp + "_" + tree_name.replace(".nhx", '.tex'), index=False, escape=False,
                    float_format=lambda x: ("{0:.3g}".format(x) if is_float(x) else x))
