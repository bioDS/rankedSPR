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
import simulate_trees as sim


def all_pw_dist(input_file, output_file = '', distances_file = '', metric = 'RNNI'):
    # plot and save histogram of all pw distances ('RNNI' od 'RF'). inpput_file is filehandle to file with input trees, distances_file (if specified) is either file to read the matrix from (if already exists), or file to save matrix in (make sure it didn't exist before!), output_file is filehandle for saving the plot (histoogram).
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation)
        print("Read trees")
        tree_list = read_nexus(input_file, ranked = True)
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
        tree_list = read_nexus(input_file, ranked = True)
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
        tree_list = read_nexus(input_file, ranked = True)
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

        # Plotting RF distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
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
        tree_list = read_nexus(input_file, ranked = True)
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RNNI distances
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

        # Plotting RF distances
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
        init_tree_list = read_nexus(input_file, ranked = True)
        index = init_tree_list.num_trees # index that the focal tree will have in tree_list (focal tree is last tree in that list)
        trees = (TREE * (index + 1))() # list of trees containing focal tree as last tree
        for i in range(0, index):
            trees[i] = init_tree_list.trees[i]
        trees[index] = focal_tree
        tree_list = TREE_LIST(trees, index + 1)
        tree_list.num_trees = index + 1
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, the last tree in tree_list
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
        tree_list.append(focal_tree) # add focal tree as last tree to tree_lists
        print("Done reading trees")
        num_leaves = len(leaf_labels)
        num_trees = len(tree_list)
        rf_diameter = int(2*(num_leaves - 1))

        # Plotting RF distances for all pairs T_index, T_i, where index belongs to the focal tree, the last tree in tree_list
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rf.rf_distance_focal(tree_list, len(tree_list)-1)[0]
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rf_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


# use own implementation of coalescent to plot RNNI distances
def coal_pw_dist(num_leaves, num_trees, mean = False, output_file = '', distances_file = ''):
    # Plotting the distances of all tree pairs T_i, T_i+1 for even i (for list of simulated trees this should give independent distances) and save plot (if filehandle given) in output_file
    # Read trees in C format (for RNNI distance computation)
    # If mean == True, returns mean and var of distances
    print("Simulate trees")
    tree_list = sim.sim_coal(num_leaves,num_trees)
    print("Done simulating trees")
    num_trees = tree_list.num_trees
    num_leaves = tree_list.trees[0].num_leaves
    rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

    # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
    if os.path.exists(distances_file):
        distances = np.loadtxt(distances_file, delimiter = ' ')
    else:
        distances = rnni.rnni_distances_tree_pairs(tree_list)[0]
        if mean == True:
            return(np.mean(distances), np.var(distances))
        if distances_file != '':
            np.savetxt(distances_file, distances, delimiter = ' ')
    if mean == False:
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


# simulate coalescent trees and plot distance to caterpillar tree [1,2,...,n]
def caterpillar_dist_distribution(num_leaves, num_trees, mean = False, output_file = '', distances_file = ''):
    # If mean == True, returns mean and var of distances
    sim_trees = sim.sim_coal(num_leaves, num_trees).trees
    all_trees = (TREE * (num_trees + 1))()
    rnni_diameter = (num_leaves - 1)*(num_leaves - 2)/2
    for i in range(0,num_trees):
        all_trees[i] = sim_trees[i]
    all_trees[num_trees] = sim.identity_caterpillar(num_leaves)
    tree_list = TREE_LIST(all_trees, num_trees)
    # Plotting RNNI distances for all pairs T_{num_trees}, T_i, where num_trees belongs to the focal tree, which is the caterpillar trees 1,2,...,n
    if os.path.exists(distances_file):
        distances = np.loadtxt(distances_file, delimiter = ' ')
    else:
        distances = rnni.rnni_distance_focal(tree_list, num_trees)[0]
        norm_distances = []
        for i in distances:
            if rnni_diameter != 0:
                norm_distances.append(i/rnni_diameter)
        if mean == True:
            return(np.mean(norm_distances), np.var(norm_distances))
        if distances_file != '':
            np.savetxt(distances_file, distances, delimiter = ' ')
    if mean == False:
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


# use own implementation of coalescent to plot RNNI distances
def coal_focal_dist(num_leaves, num_trees, mean = False, output_file = '', distances_file = ''):
    # Simulating num_trees coalescent trees and plotting the distances of all tree pairs T_i, T_{num_trees} for all i and save plot (if filehandle given) in output_file
    # If mean == True, returns mean and var of distances
    # Read trees in C format (for RNNI distance computation)
    print("Simulate trees")
    tree_list = sim.sim_coal(num_leaves,num_trees)
    print("Done simulating trees")
    num_trees = tree_list.num_trees
    num_leaves = tree_list.trees[0].num_leaves
    rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

    # Plotting RNNI distances for all pairs T_{num_trees}, T_i, where num_trees belongs to the focal tree
    if os.path.exists(distances_file):
        distances = np.loadtxt(distances_file, delimiter = ' ')
    else:
        distances = rnni.rnni_distance_focal(tree_list, num_trees-1)[0]
        norm_distances = []
        for i in distances:
            if rnni_diameter != 0:
                norm_distances.append(i/rnni_diameter)
        if mean == True:
            return(np.mean(norm_distances), np.var(norm_distances))
        if distances_file != '':
            np.savetxt(distances_file, distances, delimiter = ' ')
    if mean == False:
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


def dist_distribution_to_caterpillars(num_leaves, num_trees, mean = False, output_file = '', distances_file = ''):
    # Simulate random tree and compute distance to num_trees simulated caterpillar trees
    # If mean == True, returns mean and var of distances
    cat_trees = sim.sim_cat(num_leaves, num_trees).trees
    all_trees = (TREE * (num_trees + 1))()
    rnni_diameter = (num_leaves - 1)*(num_leaves - 2)/2
    for i in range(0,num_trees):
        all_trees[i] = cat_trees[i]
    all_trees[num_trees] = sim.sim_coal(num_leaves, 1).trees[0] # simulated coalescence tree is focal tree here
    tree_list = TREE_LIST(all_trees, num_trees)
    # Plotting RNNI distances for all pairs T_{num_trees}, T_i, where num_trees belongs to the focal tree, a tree sampled from coalescent
    if os.path.exists(distances_file):
        distances = np.loadtxt(distances_file, delimiter = ' ')
    else:
        distances = rnni.rnni_distance_focal(tree_list, num_trees)[0]
        norm_distances = []
        for i in distances:
            if rnni_diameter != 0:
                norm_distances.append(i/rnni_diameter)
        if mean == True:
            return(np.mean(norm_distances), np.var(norm_distances))
        if distances_file != '':
            np.savetxt(distances_file, distances, delimiter = ' ')
    if mean == False:
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


def dist_distribution_btw_caterpillars(num_leaves, num_trees, mean = False, output_file = '', distances_file = ''):
    # Simulate random tree and compute distance to num_trees simulated caterpillar trees
    # If mean == True, returns mean and var of distances
    trees = sim.sim_cat(num_leaves, num_trees).trees
    rnni_diameter = (num_leaves - 1)*(num_leaves - 2)/2
    tree_list = TREE_LIST(trees, num_trees)
    # Plotting RNNI distances for all pairs T_{num_trees}, T_i, where num_trees belongs to the focal tree, a tree sampled from coalescent
    if os.path.exists(distances_file):
        distances = np.loadtxt(distances_file, delimiter = ' ')
    else:
        distances = rnni.rnni_distances_tree_pairs(tree_list, num_trees)[0]
        norm_distances = []
        for i in distances:
            if rnni_diameter != 0:
                norm_distances.append(i/rnni_diameter)
        if mean == True:
            return(np.mean(norm_distances), np.var(norm_distances))
        if distances_file != '':
            np.savetxt(distances_file, distances, delimiter = ' ')
    if mean == False:
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


def mean_distance_n(func, min_num_leaves, max_num_leaves, num_trees, output_file = ''):
    # plot mean distances given by function for different number of leaves and plot them
    mean_array = []
    var_array = []
    for i in range(min_num_leaves,max_num_leaves):
        statistics = func(i,20000, mean = True)
        mean_array.append(statistics[0])
        var_array.append(statistics[1])
    print(mean_array)
    print(var_array)
    plt.plot(mean_array)
    plt.plot(var_array)
    if output_file != '':
        plt.savefig(output_file)
    plt.show()


def mean_distance_repeat(func, num_leaves, num_iterations, num_trees, output_file = ''):
    # plot mean distances given by function for different number of leaves and plot them
    mean_array = []
    var_array = []
    for i in range(0,num_iterations):
        statistics = func(num_leaves,num_trees, mean = True)
        mean_array.append(statistics[0])
        var_array.append(statistics[1])
    print(mean_array)
    print(var_array)
    plt.plot(mean_array)
    plt.plot(var_array)
    if output_file != '':
        plt.savefig(output_file)
    plt.show()


if __name__ == '__main__':

    # dist_distribution_to_caterpillars(20,10000, output_file = '../simulations/distance_distribution/coalescent/dist_distribution_to_caterpillars_20_n_10000_N.eps')
    # dist_distribution_btw_caterpillars(20,20000, output_file = '../simulations/distance_distribution/coalescent/dist_distribution_btw_caterpillars_20_n_20000_N.eps')
    mean_distance_n(dist_distribution_btw_caterpillars, 4, 40, 10000, output_file = '../simulations/distance_distribution/coalescent/btw_cat_mean_and_var_dist_n_3_to_40_N_20000.eps')
    # mean_distance_repeat(dist_distribution_to_caterpillars, 20, 50, 1000, output_file = '../simulations/distance_distribution/coalescent/to_cat_mean_and_var_dist_n_20_N_20000_50_iterations.eps')
    # coal_pw_dist(20000,20, output_file = '../simulations/distance_distribution/coalescent/own_coal_distr_20_n_20000_N.eps')
    # caterpillar_dist_distribution(20,20000, output_file='../simulations/distance_distribution/coalescent/caterpillar_distances_20_n_20000_N.eps')
    # coal_focal_dist(20, 20000, output_file='../simulations/distance_distribution/coalescent/coal_focal_dist_20_n_20000_N.eps')

    # pw_tree_list_dist('../simulations/simulated_trees/bd/20000/bd_trees_20_n.nex', '../simulations/distance_distribution/bd/rnni_distribution_20_n_20000_N.eps', metric = 'RNNI')
    # pw_tree_list_dist('../simulations/posterior/coal/coal_alignment_20_sequences_10000_length.trees', '../simulations/posterior/coal/rf_all_pw_dist.eps', metric = 'RF')
    # read MCC tree:

    # mcc_ete_tree = rf.read_ete_nexus('../simulations/posterior/bd/mcc_summary.tree')[0]
    # mcc_ete_tree = mcc_ete_tree[0]
    # focal_tree_dist(mcc_ete_tree, input_file = '../simulations/posterior/bd/bd_alignment_b1_d0_20_sequences_10000_length.trees', output_file = '../simulations/posterior/bd/rf_mcc_dist.eps', metric = 'RF')

    # mcc_tree = read_nexus('../simulations/posterior/coal/mcc_summary.tree', ranked = True).trees[0]
    # focal_tree_dist(mcc_tree, '../simulations/posterior/coal/coal_alignment_20_sequences_10000_length.trees', '../simulations/posterior/coal/rnni_mcc_dist.eps', metric = 'RNNI')

    # f = open('../simulations/posterior/coal/coal_20_leaves.new', 'r')
    # original_tree = read_newick(f.readline(), ranked = True)
    # print(findpath_distance(original_tree, mcc_tree))
    # focal_tree_dist(original_tree, '../simulations/posterior/coal/coal_alignment_20_sequences_10000_length.trees', '../simulations/posterior/coal/rnni_original_tree_dist.eps', metric = 'RNNI')
    # f.close()