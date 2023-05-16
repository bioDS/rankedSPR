__author__ = 'Lena Collienne'
# Functions used for ranked SPR paper
import sys

from rankedSPR_seidel import *
from rankedspr_exploration import *
from spr_path_functions import *


def get_diameter(num_leaves):
    d = get_distance_matrix(num_leaves, hspr=True)
    print("The diameter of HSPR space for", num_leaves, "leaves is", np.amax(d[0]))

def check_approx_alg(num_leaves):
    (d, tree_dict, tree_index_dict) = get_distance_matrix(num_leaves,
                                                           hspr=True)
    for i in range(0, len(d)):
        tree1 = read_from_cluster(tree_index_dict[i])
        for j in range(i, len(d)):
            tree2 = read_from_cluster(tree_index_dict[j])
            approx_path = rankedspr_path_bottom_up_hspr(tree1, tree2)
            approx_dist = approx_path.num_trees - 1
            if approx_dist != d[i][j]:
                print("For the trees", tree_index_dict[i], "and", tree_index_dict[j],
                    "the approximated HSPR distance is", approx_dist)
                print("The exact distance is", d[i][j])
                return(False)
    return(True)

def main():
    n = 4
    get_diameter(n)
    check_approx_alg(n)

    tree1_str = "[{4,5}:1,{1,2}:2,{1,2,3}:3,{1,2,3,4,5}:4]"
    tree2_str = "[{1,2}:1,{4,5}:2,{1,2,3}:3,{1,2,3,4,5}:4]"
    tree1 = read_from_cluster(tree1_str)
    tree2 = read_from_cluster(tree2_str)

    # dictionary containing predecessors of tree2 for all paths from tree1 to tree2
    pred_dict = all_shortest_paths(tree1, tree2)

    # check if any tree in the predecessor dictionary (other than tree1 and tree2)
    # contains the cluster {1,2,3}
    cluster_found = False
    for tree_str in pred_dict:
        if tree_str != tree1_str and tree_str != tree2_str and "{1,2,3}" in tree_str:
            print(tree_str + " contains the cluster {1,2,3}")
            cluster_found = True
        for pred_str in pred_dict[tree_str]:
            if pred_str != tree1_str and pred_str != tree2_str and "{1,2,3}" in pred_str:
                print(tree_str + " contains the cluster {1,2,3}")
                cluster_found = True
        if cluster_found:
            break
    if cluster_found:
        print("The shared clusters {1,2,3} is present in a tree on a shortes path between "
        + tree1_str + " and " + tree2_str)
    else:
        print("The shared clusters {1,2,3} is not present in a tree on a shortes path between "
        + tree1_str + " and " + tree2_str)


if __name__ == "__main__":
    main()
