__author__ = 'Lena Collienne'

import os.path
import sys
sys.path.append('../..')

from ete3 import Tree
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from tree_structs import *
from dct_parser.tree_io import *
import rf_distances as rf
import rnni_distances as rnni
import plots as plts


def all_pw_dist(input_file, output_file = '', distances_file = '', metric = 'RNNI'):
    # plot and save histogram of all pw distances ('RNNI' od 'RF'). inpput_file is filehandle to file with input trees, distances_file (if specified) is either file to read the matrix from (if already exists), or file to save matrix in (make sure it didn't exist before!), output_file is filehandle for saving the plot (histoogram).
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation)
        print("Read trees")
        tree_list = read_nexus(input_file, ranked = True)[0]
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RNNI distances for all pairs T_i, T_j, i<j
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rnni.pw_rnni_dist(tree_list)
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)
    
    elif metric == 'RF':
        # Read trees in ete3 format:
        print("Read trees")
        tree_list, leaf_labels = rf.read_ete_nexus(input_file)
        print("Done reading trees")
        num_leaves = len(leaf_labels)
        rf_diameter = int(2*(num_leaves - 1))

        # Plotting RF distances for all pairs T_i, T_j, i<j
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rf.pw_rf_dist(tree_list, list = True)
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rf_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)

def random_focal_tree_dist(input_file, output_file = '', distances_file = '', metric = 'RNNI'):
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation)
        print("Read trees")
        tree_list = read_nexus(input_file, ranked = True)[0]
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            index = np.random.randint(0,num_trees)
            distances = rnni.rnni_distance_focal(tree_list, index)[0]
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)

    elif metric == 'RF':
        # Read trees in ete3 format:
        print("Read trees")
        tree_list, leaf_labels = rf.read_ete_nexus(input_file)
        print("Done reading trees")
        num_leaves = len(leaf_labels)
        num_trees = len(tree_list)
        rf_diameter = int(2*(num_leaves - 1))

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            index = np.random.randint(0,num_trees)
            distances = rf.rf_distance_focal(tree_list, index)[0]
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rf_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


def consec_trees_dist(input_file, output_file = '', distances_file = '', metric = 'RNNI'):
    # Plotting the distances of all tree pairs T_i, T_i+1 and save plot (if filehandle given) in output_file
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation)
        print("Read trees")
        tree_list = read_nexus(input_file, ranked = True)[0]
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rnni.rnni_distances_consecutive_tree_pairs(tree_list)[0]
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        plts.plot_dots(distances, output_file)

    elif metric == 'RF':
        # Read trees in ete3 format:
        print("Read trees")
        tree_list, leaf_labels = rf.read_ete_nexus(input_file)
        print("Done reading trees")
        num_leaves = len(leaf_labels)
        num_trees = len(tree_list)
        rf_diameter = int(2*(num_leaves - 1))

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rf.rf_distances_consecutive_tree_pairs(tree_list)[0]
            print(distances)
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        plts.plot_dots(distances, output_file)


def pw_tree_list_dist(input_file, output_file = '', distances_file = '', metric = 'RNNI'):
    # Plotting the distances of all tree pairs T_i, T_i+1 for even i (for list of simulated trees this should give independent distances) and save plot (if filehandle given) in output_file
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation)
        print("Read trees")
        tree_list = read_nexus(input_file, ranked = True)[0]
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rnni.rnni_distances_tree_pairs(tree_list)[0]
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)

    elif metric == 'RF':
        # Read trees in ete3 format:
        print("Read trees")
        tree_list, leaf_labels = rf.read_ete_nexus(input_file)
        print("Done reading trees")
        num_leaves = len(leaf_labels)
        num_trees = len(tree_list)
        rf_diameter = int(2*(num_leaves - 1))

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rf.rf_distances_tree_pairs(tree_list)[0]
            print(distances)
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rf_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)

def focal_tree_dist(focal_tree, input_file, output_file = '', distances_file = '', metric = 'RNNI'):
    # Compute distances from focal_tree to all trees in input_file and save as histogram
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation) -- does not contain focal tree yet.
        print("Read trees")
        init_tree_list = read_nexus(input_file, ranked = True)[0]
        index = init_tree_list.num_trees # index that the focal tree will have in tree_list (focal tree is last tree in that list)
        trees = (TREE * (index + 1))() # list of trees containing focal tree as last tree
        for i in range(0, index):
            trees[i] = init_tree_list.trees[i]
        trees[index] = focal_tree
        tree_list = TREE_LIST(trees, index + 1)
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rnni.rnni_distance_focal(tree_list, index)[0]
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)

    elif metric == 'RF':
        # Read trees in ete3 format:
        print("Read trees")
        tree_list, leaf_labels = rf.read_ete_nexus(input_file)
        # Correct the leaf labels in focal tree:
        inv_leaf_labels = dict() # inverse of leaf_label dict
        for key in leaf_labels.keys():
            inv_leaf_labels[leaf_labels[key]] = key
        for node in mcc_ete_tree.get_leaves(): # Replace leaf labels
            if "t" in node.name:
                newname = inv_leaf_labels[node.name]
                node.name = newname
        tree_list.append(focal_tree) # add focal tree as last tree to tree_lists
        print("Done reading trees")
        num_leaves = len(leaf_labels)
        num_trees = len(tree_list)
        rf_diameter = int(2*(num_leaves - 1))

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rf.rf_distance_focal(tree_list, len(tree_list)-1)[0]
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rf_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


if __name__ == '__main__':

    # all_pw_dist('../simulations/posterior/coal/coal_alignment_20_sequences_10000_length.trees', '../simulations/posterior/coal/rnni_all_pw_dist.eps', metric = 'RNNI')
    # all_pw_dist('../simulations/posterior/coal/coal_alignment_20_sequences_10000_length.trees', '../simulations/posterior/coal/rf_all_pw_dist.eps', metric = 'RF')
    # read MCC tree:
    f = open('../simulations/posterior/coal/mcc_summary.new')
    tree_str = f.readline()
    # mcc_tree = read_newick(tree_str, ranked = True)
    mcc_ete_tree = Tree(tree_str)
    # focal_tree_dist(mcc_tree, '../simulations/posterior/coal/coal_alignment_20_sequences_10000_length.trees', '../simulations/posterior/coal/rnni_mcc_dist.eps', metric = 'RNNI')
    focal_tree_dist(mcc_ete_tree, input_file = '../simulations/posterior/coal/coal_alignment_20_sequences_10000_length.trees', output_file = '../simulations/posterior/coal/rf_mcc_dist.eps', metric = 'RF')
    
    f.close()