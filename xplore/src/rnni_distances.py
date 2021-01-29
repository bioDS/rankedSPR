__author__ = 'Lena Collienne'

import sys
sys.path.append('../..')

from ete3 import Tree
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from tree_structs import *
from dct_parser.tree_io import *
import plots as plts


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


def rnni_distances_consecutive_tree_pairs(tree_list, prgr = False):
    # returns array of distances between every pair of trees i,i+1 in given file and number of leaves
    distances = []
    if (prgr == True):
        progress = 0.05 #for printing progress
    print("Computing RNNI distances")
    for i in range(0, tree_list.num_trees - 1):
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


def pw_rnni_dist(tree_list):
    # Return either a np.matrix of pw distances (if list=False; this is an upper diagonal matrix!), or a sequence containing distances (excluding diagonal entries)
    num_trees = tree_list.num_trees
    if list == False:
        # Create empty distance matrix as input for pw_distances
        distances = np.zeros(shape=(num_trees,num_trees),dtype=np.int32)
    else:
        distances = []
    for i in range(0,num_trees):
        for j in range(i + 1,num_trees):
            if list == False:
                distances[i][j] = (findpath_distance(tree_list.trees[i],tree_list.trees[j]))
            else:
                distances.append(findpath_distance(tree_list.trees[i],tree_list.trees[j]))
    return distances


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


# if __name__ == '__main__':

    # plt.plot(rnni_mean_dist_n(100, 20000), linestyle = 'None', marker = 'o', markersize = 6) # coalescent
    # plt.plot(rnni_mean_dist_n(40, 20000, model = 'bd'), linestyle = 'None', marker = 'o', markersize = 6) # birth-death
    # plt.show()

    # Get input trees:
    # filename = input("What is the file with trees?\n")

    # # Read trees in C format (for RNNI distance computation)
    # print("Read trees")
    # tree_list = read_nexus(filename, ranked = True)[0]
    # print("Done reading trees")
    # num_trees = tree_list.num_trees
    # num_leaves = tree_list.trees[0].num_leaves
    # rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

    # # Plotting RNNI (consecutive pairs) distances
    # print(num_leaves)
    # distances_rnni,num_leaves = rnni_distances_consecutive_tree_pairs(tree_list)
    # print(np.mean(distances_rnni))
    # print(distances_rnni)
    # plt.plot(distances_rnni)
    # plt.show()

    # # all PW RNNI distances:
    # pw_distances_rnni = pw_rnni_dist(tree_list)
    # plt.hist(pw_distances_rnni, bins = rnni_diameter, range = (0, rnni_diameter))
    # plt.show()

    # # Plotting RNNI pw distances (between pairs of trees T_i, T_{i+1} for even i)
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
    # tree_list, leaf_labels = rf.read_ete_nexus(filename)
    # num_leaves = len(leaf_labels)

    # # all PW RF distances:
    # pw_distances_rf = rf.pw_rf_dist(tree_list, list = True)

    # # Save/load pw distance matrix:
    # np.savetxt('../simulations/posterior/coal/pw_RF_distance_matrix.txt', pw_distances_rf, delimiter = ' ')
    # pw_distances_rf = np.loadtxt('../simulations/posterior/coal/pw_RF_distance_matrix.txt', delimiter = ' ')
    # plts.plot_hist(pw_distances_rf, '../simulations/posterior/coal/rf_all_pw_dist.eps')

    # # Plotting RF distances:
    # distances_rf, num_leaves = rf.rf_distances_tree_pairs(tree_list)
    # rf_diameter = int(2*(num_leaves - 1))
    # plt.hist(distances_rf, bins = rf_diameter, range = (0, rf_diameter))
    # plt.show()

    # # Plotting RF (consecutive pairs) distances:
    # distances_rf, num_leaves = rf.rf_distances_consecutive_tree_pairs(tree_list)
    # rf_diameter = int(2*(num_leaves - 1))
    # plt.plot(distances_rf)
    # plt.show()
    # plt.hist(distances_rf, bins = rf_diameter, range = (0, rf_diameter))
    # plt.show()

    # # Plotting RNNI distances to random focal tree
    # num_trees = len(tree_list)
    # index = np.random.randint(0,num_trees)
    # focal_dist, num_leaves = rf.rf_distance_focal(tree_list, index)
    # rf_diameter = int(2*(num_leaves - 1))
    # plt.hist(focal_dist, bins = rf_diameter ,range = (0, rf_diameter))
    # plt.show()
