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
        plts.plot_hist(distances, output_file, bins = rnni_diameter)
    
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
        plts.plot_hist(distances, output_file, bins = rf_diameter)

def focal_tree_dist(input_file, output_file = '', distances_file = '', metric = 'RNNI'):
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
        plts.plot_hist(distances, output_file, bins = rnni_diameter)

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
        plts.plot_hist(distances, output_file, bins = rf_diameter)


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
        plts.plot_dots(distances, output_file)


if __name__ == '__main__':

    pw_tree_list_dist('../simulations/simulated_trees/coal/20000/coal_trees_6_n.nex', '../simulations/distance_distribution/coalescent/rnni_distribution_6_n_20000_N.eps', metric = 'RNNI')
    # consec_trees_dist('../simulations/simulated_trees/coal/coal_trees_20_n_100_N.nex', '../simulations/simulated_trees/coal/output.eps', metric = 'RNNI')
