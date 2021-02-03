__author__ = 'Lena Collienne'

import sys
sys.path.append('../..')
from tree_structs import *
import re
from ete3 import Tree
import matplotlib.pyplot as plt
import numpy as np


def read_ete_nexus(file_handle):
    # Read trees from nexus file and save as list fo ete tree
    # To get the number of trees in the given file, we find the first and last line in the file that contain a tree (Line starts with tree)
    # Count number of lines in file
    last_line = len(open(file_handle).readlines())

    # Find last line containing a tree
    for line in reversed(list(open(file_handle))):
        re_tree = re.search(r'tree', line, re.I)
        if re_tree == None:
            last_line -= 1
        else: break

    # Find first line containing a tree
    first_line = 1
    in_tree = False # turn to True if we pass 'begin tree' in NEXUS file. The actual trees then start after the 'translate' block
    for line in list(open(file_handle)):
        if in_tree == False:
            re_tree = re.search(r'begin tree', line, re.I)
            if re_tree != None:
                in_tree = True
            first_line += 1
        else: # pass the translate block
            if re.search(r';', line) == None:
                first_line += 1
            else:
                first_line += 1
                break

    num_trees = last_line - first_line + 1 # Number of trees in nexus file

    # running variables for reading trees and displaying progress
    index = 0
    progress = 10

    name_dict = dict() # Save tree label names in dict
    trees = list()

    f = open(file_handle, 'r')

    leaf_labels = False
    # Save leaf labels -- returns empty dict if no abbreviation for leaf labels is used
    for line in f:
        if re.search(r'translate', line, re.I) != None:
            leaf_labels = True
        # Start reading leaf labels after 'translate' (signals start of this sequence)
        if leaf_labels == True:
            re_label = re.search(r'\s*(.+)\s([^,^\n]+)', line)
            re_stop = re.search(r';', line)
            if re_stop != None:
                break
            elif re_label != None:
                name_dict[re_label.group(1)] = re_label.group(2)

    # Read trees
    for line in f:
        re_tree = re.search(r'tree .* (\(.*\);)', line, re.I)
        if re_tree != None:
            current_tree = Tree(re.sub(r'\[[^\]]*\]',"",re_tree.group(1)))
            trees.append(current_tree)          
            index += 1
            if int(100*index/num_trees) == progress:
                print(str(progress) + '% of trees are read')
                progress += 10
    f.close()
    return(trees, name_dict)


def rf_distances_tree_pairs(tree_list):
    # returns array of distances between every pair of trees i,i+1 (even i) in given file and number of leaves
    num_trees = len(tree_list)
    print('number of trees: ', num_trees)
    distances = []
    progress = 0.05 #for printing progress
    num_leaves = len(tree_list[0])
    print("Computing RF distances")
    for i in range(0, num_trees - 1, 2):
        distances.append(tree_list[i].robinson_foulds(tree_list[i+1])[0])
        if (i/(num_trees) > progress):
            print('Progress: ' + "{:.2f}".format(progress))
            progress += 0.05
    return(distances, num_leaves)


def rf_distances_consecutive_tree_pairs(tree_list):
    # returns array of distances between every pair of trees i,i+1 in given file and number of leaves
    num_trees = len(tree_list)
    print('number of trees: ', num_trees)
    distances = []
    progress = 0.05 #for printing progress
    num_leaves = len(tree_list[0])
    print("Computing RF distances")
    for i in range(0, num_trees - 1):
        distances.append(tree_list[i].robinson_foulds(tree_list[i+1])[0])
        if (i/(num_trees) > progress):
            print('Progress: ' + "{:.2f}".format(progress))
            progress += 0.05
    return(distances, num_leaves)


def rf_distance_focal(tree_list, index):
    # returns array of distances from chosen focal tree (tree number i in nexus file)
    distances = []
    progress = 0.05 #for printing progress
    num_trees = len(tree_list)
    num_leaves = len(tree_list[0])
    print("Computing RNNI distances")
    for i in range(0, num_trees):
        if (i != index):
            distances.append(tree_list[i].robinson_foulds(tree_list[index])[0])
        if (i/num_trees > progress):
            print('Progress: ' + "{:.2f}".format(progress))
            progress += 0.05
    return(distances, int(num_leaves))


def pw_rf_dist(tree_list, list = False):
    # Return either a np.matrix of pw distances (if list=False; this is an upper diagonal matrix!), or a sequence containing distances (excluding diagonal entries)
    num_trees = len(tree_list)
    if list == False:
        distances = np.zeros(shape=(num_trees,num_trees),dtype=np.int32)
    else:
        distances = []
    for i in range(0,num_trees):
        for j in range(i + 1,num_trees):
            if list == False:
                distances[i][j] = (tree_list[i].robinson_foulds(tree_list[j])[0])
            else:
                distances.append(tree_list[i].robinson_foulds(tree_list[j])[0])
    return distances