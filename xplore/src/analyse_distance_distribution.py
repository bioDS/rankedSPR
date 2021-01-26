__author__ = 'Lena Collienne'

import sys
sys.path.append('../..')

from ete3 import Tree
import matplotlib.pyplot as plt
import numpy as np

from tree_structs import *
from dct_parser.tree_io import *
import rf_distance_analysis as rf


def rnni_distances_tree_pairs(tree_list, prgr = False):
    # returns array of distances between every pair of trees i,i+1 (even i) in given file and number of leaves
    distances = []
    if (prgr == True):
        progress = 0.05 #for printing progress
    print("Computing RNNI distances")
    for i in range(0, tree_list.num_trees, 2):
        # print(i,i+1)
        distances.append(findpath_distance(tree_list.trees[i],tree_list.trees[i+1]))
        if prgr == True and (i/(tree_list.num_trees) > progress):
            print('Progress: ' + "{:.2f}".format(progress))
            progress += 0.05
        # print(i, tree_to_cluster_string(tree_list.trees[i]), i+1, tree_to_cluster_string(tree_list.trees[i+1]))
    return(distances, int(tree_list.trees[0].num_leaves))

def rnni_distance_focal(tree_list, index, prgr = False):
    # returns array of distances from chosen focal tree (tree number i in nexus file)
    distances = []
    if (prgr == True):
        progress = 0.05 #for printing progress
    print("Computing RNNI distances")
    for i in range(0, tree_list.num_trees):
        if (i != index):
            distances.append(findpath_distance(tree_list.trees[index],tree_list.trees[i]))
        if prgr == True and (i/(tree_list.num_trees) > progress):
            print('Progress: ' + "{:.2f}".format(progress))
            progress += 0.05
    return(distances, int(tree_list.trees[0].num_leaves))

def rnni_mean_dist_n(n,N,model='coal'):
    # plots the mean distance between trees given in a file in the folder (as computed in rnni_distances_tree_pairs) against n (for varying n)
    mean_list = []
    for i in range(3,n):
        print('Trees on '+ str(i) + ' leaves')
        if model == 'coal':
            filename = '../simulated_trees/coal/' + str(N) + '/coal_trees_' + str(i) + '_n.nex'
        elif model == 'bd':
            filename = '../simulated_trees/bd/' + str(N) + '/bd_trees_' + str(i) + '_n.nex'
        tree_list = read_nexus(filename, ranked = True)[0]
        distances, num_leaves = rnni_distances_tree_pairs(tree_list)
        diameter = (num_leaves-1)*(num_leaves-2)/2
        mean_list.append(np.mean(distances)/diameter)
    return(mean_list)


if __name__ == '__main__':

    # plt.plot(rnni_mean_dist_n(40, 20000), linestyle = 'None', marker = 'o', markersize = 6) # coalescent
    plt.plot(rnni_mean_dist_n(40, 20000, model = 'bd'), linestyle = 'None', marker = 'o', markersize = 6) # birth-death
    plt.show()

    # Get input trees:
    # filename = input("What is the file with trees?\n")
    # # Read trees in C format (for RNNI distance computation)
    # print("Read trees")
    # tree_list = read_nexus(filename, ranked = True)[0]
    # print("Done reading trees")
    # num_trees = tree_list.num_trees
    # num_leaves = tree_list.trees[0].num_leaves
    # rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

    # Plotting RNNI distances
    # distances_rnni,num_leaves = rnni_distances_tree_pairs(tree_list)
    # print(np.mean(distances_rnni))
    # plt.hist(distances_rnni, bins = rnni_diameter, range = (0, rnni_diameter))
    # plt.show()

    # # Plotting RNNI distances to random focal tree
    # index = np.random.randint(0,num_trees)
    # focal_dist = rnni_distance_focal(tree_list, index)
    # plt.hist(focal_dist, bins = rnni_diameter ,range = (0, rnni_diameter))
    # plt.show()

    # # Read trees in ete3 format (for RF distance)
    # tree_list, num_leaves = rf.read_ete_nexus(filename)

    # # Plotting RF distances:
    # distances_rf, num_leaves = rf.rf_distances_tree_pairs(filename)
    # rf_diameter = int(2*(num_leaves - 1))
    # plt.hist(distances_rf, bins = rf_diameter, range = (0, rf_diameter))
    # plt.show()

    # # Plotting RNNI distances to random focal tree
    # num_trees = len(tree_list)
    # index = np.random.randint(0,num_trees)
    # focal_dist, num_leaves = rf.rf_distance_focal(tree_list, index)
    # rf_diameter = int(2*(num_leaves - 1))
    # plt.hist(focal_dist, bins = rf_diameter ,range = (0, rf_diameter))
    # plt.show()