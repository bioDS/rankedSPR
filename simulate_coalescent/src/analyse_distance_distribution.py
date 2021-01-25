__author__ = 'Lena Collienne'

import sys
sys.path.append('../..')

from ete3 import Tree
import matplotlib.pyplot as plt
import numpy as np

from tree_structs import *
from dct_parser.tree_io import *
import rf_distance_analysis as rf


def rnni_distances_tree_pairs(filename):
    # returns array of distances between every pair of trees i,i+1 (even i) in given file and number of leaves
    print("Read trees")
    tree_list = read_nexus(filename, ranked = True)[0]
    distances = []
    progress = 0.05 #for printing progress
    print("Computing RNNI distances")
    for i in range(0, tree_list.num_trees, 2):
        # print(i,i+1)
        distances.append(findpath_distance(tree_list.trees[i],tree_list.trees[i+1]))
        if (i/(tree_list.num_trees) > progress):
            print('Progress: ' + "{:.2f}".format(progress))
            progress += 0.05
        # print(i, tree_to_cluster_string(tree_list.trees[i]), i+1, tree_to_cluster_string(tree_list.trees[i+1]))
    return(distances, int(tree_list.trees[0].num_leaves))

def rnni_distance_focal(filename, index):
    # returns array of distances from chosen focal tree (tree number i in nexus file)
    print("Read trees")
    tree_list = read_nexus(filename, ranked = True)[0]
    distances = []
    progress = 0.05 #for printing progress
    print("Computing RNNI distances")
    for i in range(0, tree_list.num_trees):
        if (i != index):
            distances.append(findpath_distance(tree_list.trees[index],tree_list.trees[i]))
        if (i/(tree_list.num_trees) > progress):
            print('Progress: ' + "{:.2f}".format(progress))
            progress += 0.05
    return(distances, int(tree_list.trees[0].num_leaves))


if __name__ == '__main__':
    # Get input trees:
    tree_file = input("What is the file with trees?\n")
    # # Plotting RNNI distances
    # distances_rnni,num_leaves = rnni_distances_tree_pairs(tree_file)
    # rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)
    # plt.hist(distances_rnni, bins = rnni_diameter, range = (0, rnni_diameter))
    # plt.show()
    # # Plot RNNI distances to random focal tree
    # num_trees = 2*len(distances_rnni)
    num_trees = 20000
    index = np.random.randint(0,num_trees)
    focal_dist = rnni_distance_focal(tree_file, index)
    plt.hist(focal_dist, bins = rnni_diameter, range = (0, rnni_diameter))
    plt.show()
    # # Plotting RF distances:
    # distances_rf, num_leaves = rf.rf_distances_tree_pairs(tree_file)
    # rf_diameter = int(2*(num_leaves - 1))
    # plt.hist(distances_rf, bins = rf_diameter, range = (0, rf_diameter))
    # plt.show()
