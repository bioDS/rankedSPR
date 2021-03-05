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
    # print("Computing RNNI distances")
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
    # print("Computing RNNI distances")
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
    # print("Computing RNNI distances")
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