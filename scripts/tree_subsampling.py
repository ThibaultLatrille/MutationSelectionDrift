from ete3 import Tree
import matplotlib.pyplot as plt

current_dir = "/home/thibault/SimuEvol/"
folder = "data_trees"
protein = "np"

tree_path = "{0}/{1}/{2}.newick".format(current_dir, folder, protein)

t = Tree(tree_path)
min_n = 25
nbr_leaves, distance_farthest = [], []
for index, node in enumerate(t.traverse("levelorder")):
    n = len(node.get_leaves())
    if n > min_n:
        print(index)
        nbr_leaves.append(n)
        distance_farthest.append(node.get_farthest_leaf()[1])

leaves_range = range(min_n, max(nbr_leaves) + 1)
nbr_subtree = []
for n in leaves_range:
    nbr_subtree.append(len([nbr for nbr in nbr_leaves if nbr <= n]))

plt.plot(leaves_range, nbr_subtree)
plt.xlabel("Number of leaves")
plt.ylabel("Number of subtrees")
plt.title("Number of leaves for each subtree")
plt.show()

plt.scatter(nbr_leaves, distance_farthest)
plt.xlabel("Number of leaves")
plt.ylabel("Distance to farthest leaf")
plt.title("Distance to farthest leaf as a function of the number of leaves")
plt.show()

print('Job completed')
