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
                print("For the trees", tree_index_dict[i], "and", tree_index_dict[j], "the approximated HSPR distance is", approx_dist)
                print("The exact distance is", d[i][j])
                return(False)
    return(True)

def main():
    n = 4
    get_diameter(n)
    check_approx_alg(n)

if __name__ == "__main__":
    main()
