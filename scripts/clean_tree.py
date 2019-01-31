import os
from ete3 import Tree

current_dir = "/home/thibault/SimuEvol/"
folder = "data_primates"
folder_path = "{0}/{1}".format(current_dir, folder)
tree_extension = ".rootedtree"
tree_files = [f for f in os.listdir(folder_path) if tree_extension == f.strip()[-len(tree_extension):]]
print("Found {0} tree files".format(len(tree_files)))
hit = 0

for tree_file in sorted(tree_files):
    hit += 1
    tree_path = "{0}/{1}".format(folder_path, tree_file)
    newick_path = "{0}/{1}".format(folder_path, tree_file.replace(tree_extension, ".newick"))

    t = Tree(tree_path)
    name = 0
    for node in t.traverse():
        if node.name == "":
            name += 1
            node.name = "Node" + str(name)
        print(node.name)
    t.write(format=4, outfile=newick_path)
print("{0} tree files".format(hit))
print('Job completed')
