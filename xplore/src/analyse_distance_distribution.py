__author__ = 'Lena Collienne'

import os.path
from re import M
import sys
sys.path.append('../..')

from ete3 import Tree
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from fractions import Fraction


from tree_structs import *
from dct_parser.tree_io import *
import rf_distances as rf
import rnni_distances as rnni
import plots as plts
import simulate_trees as sim


def read_newick_tree_file(input):
    index = 0 #iterate through tree file
    f = open(input, 'r')
    num_trees = len(f.readlines())
    f.close()
    f = open(input, 'r')
    trees = (TREE * num_trees)()
    # Read trees in file line by line (assuming that file will only contain newick trees)
    for line in f:
        t = read_newick(line, ranked = True)
        trees[index] = t
        index += 1
    tree_list = TREE_LIST(trees, num_trees)
    f.close()
    return(tree_list)

def all_pw_dist(input, tl = False, output_file = '', distances_file = '', metric = 'RNNI'):
    # plot and save histogram of all pw distances ('RNNI' od 'RF'). input_file is filehandle to file with input trees, distances_file (if specified) is either file to read the matrix from (if already exists), or file to save matrix in (make sure it didn't exist before!), output_file is filehandle for saving the plot (histogram).
    # if tl=False, we assume input to be a filename, if tl = True, we assume input to be a TREE_LIST
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation)
        print("Read trees")
        if tl == False:
            tree_list = read_nexus(input, ranked = True)
        else:
            tree_list = input
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RNNI distances for all pairs T_i, T_j, i<j
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rnni.pw_rnni_dist(tree_list, list = True)
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        print('Done computing distance matrix.')
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)
    
    elif metric == 'RF':
        # Read trees in ete3 format:
        print("Read trees")
        if tl == False:
            tree_list, leaf_labels = rf.read_ete_nexus(input_file)
        else:
            tree_list = input
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
        d = pd.DataFrame(data =  distances)
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
        d = pd.DataFrame(data =  distances)
        plts.plot_dots(distances, output_file)


def pw_tree_list_dist(input_file, output_file = '', distances_file = '', metric = 'RNNI', nexus = True):
    # Plotting the distances of all tree pairs T_i, T_i+1 for even i (for list of simulated trees this should give independent distances) and save plot (if filehandle given) in output_file
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation)
        print("Read trees")
        if nexus == False:
            tree_list = read_newick_tree_file(input_file)
        else:
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
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        bins = np.arange(-.5, rf_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)

    elif metric == 'wRF':
        # Read trees in dendropy format:
        # print("Read trees")
        # tree_list, leaf_labels = rf.read_ete_nexus(input_file, ete = False)
        # print("Done reading trees")
        # num_leaves = len(leaf_labels)
        # num_trees = len(tree_list)
        rf_diameter = int(2*(20 - 1))

        # Plotting RF distances
        if os.path.exists(distances_file):
            distances = np.loadtxt(distances_file, delimiter = ' ')
        else:
            distances = rf.wrf_distances_tree_pairs(input_file)[0]
            if distances_file != '':
                np.savetxt(distances_file, distances, delimiter = ' ')
        # bins = np.arange(-.5, rf_diameter + 1.5, 1)
        max_dist = max(distances)
        bin_width = (max_dist - min(distances))/100 #We get 100 bars in histogram
        bins = np.arange(-.5, max_dist + 1.5, bin_width)
        d = pd.DataFrame(distances)
        plts.plot_hist(d, bins,  filehandle = output_file)


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
    # print("Simulate trees")
    tree_list = sim.sim_coal(num_leaves,num_trees)
    # print("Done simulating trees")
    num_trees = tree_list.num_trees
    num_leaves = tree_list.trees[0].num_leaves
    rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)
    # Plotting RNNI distances for all pairs T_index, T_i, where index belongs to the focal tree, which is chosen randomly (at uniform)
    if os.path.exists(distances_file):
        distances = np.loadtxt(distances_file, delimiter = ' ')
    else:
        distances = rnni.rnni_distances_tree_pairs(tree_list)[0]
        norm_distances = [i/rnni_diameter for i in distances]
        if mean == True:
            return(np.mean(norm_distances), np.var(norm_distances))
        if distances_file != '':
            np.savetxt(distances_file, distances, delimiter = ' ')
    if mean == False:
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        plts.plot_hist(distances, bins, output_file)


# use own implementation of coalescent to plot RNNI distances
def coal_pw_dist_space_efficient(num_leaves, num_tree_pairs, mean = False, output_file = '', distances_file = ''):
    # Plotting the distances of all tree pairs T_i, T_i+1 for even i (for list of simulated trees this should give independent distances) and save plot (if filehandle given) in output_file
    # Read trees in C format (for RNNI distance computation)
    # If mean == True, returns mean and var of distances
    rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)
    distances = []

    for i in range(0,int(num_tree_pairs)):
        tree_list = sim.sim_coal(num_leaves,2) # Simulate a pair of trees instead of a list with num_tree trees
        distances.append(rnni.rnni_distances_tree_pairs(tree_list)[0][0])
        norm_distances = [i/rnni_diameter for i in distances]
    if mean == True:
        return(np.mean(norm_distances), np.var(norm_distances))
    if distances_file != '':
        np.savetxt(distances_file,  distances, delimiter = ' ')
    if mean == False:
        # Plot histogram
        d = pd.DataFrame(data=distances)
        b = np.arange(-.5, rnni_diameter + 1.5, 1)
        sns.histplot(d, Color = '#b02538', Edgecolor = 'black', alpha=1, binwidth=1, binrange = [-.5,rnni_diameter+1.5], stat = 'density', legend = False)
        plt.xlabel("Distance")
        # plt.ylabel("Frequency")
        plt.ylabel("")
        plt.savefig("../simulations/distance_distribution/coalescent/rnni_distribution_" + str(num_leaves) + "_n_" + str(num_tree_pairs) + ".eps")
        plt.show()
        # plts.plot_hist(distances, bins, output_file)


# simulate coalescent trees and plot distance to a given focal tree
def given_focal_tree_dist(num_leaves, num_trees, focal_tree = None, mean = False, output_file = '', distances_file = ''):
    # If mean == True, returns mean and var of distances
    sim_trees = sim.sim_coal(num_leaves, num_trees).trees
    all_trees = (TREE * (num_trees + 1))()
    rnni_diameter = (num_leaves - 1)*(num_leaves - 2)/2
    for i in range(0,num_trees):
        all_trees[i] = sim_trees[i]
    if focal_tree == None:
        all_trees[num_trees] = sim.sim_coal(num_leaves,1).trees[0]
    else:
        all_trees[num_trees] = focal_tree
    tree_list = TREE_LIST(all_trees, num_trees)
    # Plotting RNNI distances for all pairs T_{num_trees}, T_i, where num_trees belongs to the given focal tree
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
def count_equal_trees(tree_list, output_file = '', data_file = ''):
    # Count the different number of trees in tree_list, and plot hist (number of occurrences of each tree on x axis)
    # This is to analyse the tree distribution -- if we have uniform distribution, each tree appears num_trees
    num_trees = tree_list.num_trees
    num_leaves = tree_list.trees[0].num_leaves
    distances = rnni.pw_rnni_dist(tree_list)
    values, counts = np.unique(distances, axis = 0, return_counts=True)
    d = pd.DataFrame(counts)
    # bins = np.arange(-.5, num_trees + 1.5, 1)
    plts.plot_dots(d)

# simulate coalescent trees and plot distance to a given focal tree
def compare_given_focal_tree_dist(num_leaves, num_trees, focal_tree1, focal_tree2, mean = False, output_file = ''):
    # If mean == True, returns mean and var of distances
    sim_trees = sim.sim_coal(num_leaves, num_trees).trees
    all_trees1 = (TREE * (num_trees + 1))()
    all_trees2 = (TREE * (num_trees + 1))()
    rnni_diameter = (num_leaves - 1)*(num_leaves - 2)/2
    for i in range(0,num_trees):
        all_trees1[i] = sim_trees[i]
        all_trees2[i] = sim_trees[i]
    all_trees1[num_trees] = focal_tree1
    all_trees2[num_trees] = focal_tree2
    tree_list1 = TREE_LIST(all_trees1, num_trees+1)
    tree_list2 = TREE_LIST(all_trees2, num_trees+1)
    # Plotting RNNI distances for all pairs T_{num_trees}, T_i, where num_trees belongs to the given focal tree
    distances1 = rnni.rnni_distance_focal(tree_list1, num_trees)[0]
    distances2 = rnni.rnni_distance_focal(tree_list2, num_trees)[0]
    norm_distances1 = []
    norm_distances2 = []
    for i in range(0,len(distances1)):
        if rnni_diameter != 0:
            norm_distances1.append(distances1[i]/rnni_diameter)
            norm_distances2.append(distances2[i]/rnni_diameter)
    if mean == True:
        return(np.mean(norm_distances1), np.var(norm_distances1), np.mean(norm_distances2), np.var(norm_distances2))
    if mean == False:
        bins = np.arange(-.5, rnni_diameter + 1.5, 1)
        df = pd.DataFrame(data = list(zip(distances1, distances2)), columns = ["fully balanced", "caterpillar"])
        sns.histplot(data=df, bins = bins, stat = 'density', legend = True)
        if output_file != '':
            plt.savefig(output_file)
        plt.tight_layout()
        plt.show()

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


def mean_distance_n(func, min_num_leaves, max_num_leaves, num_tree_pairs, output_file = '', variance = False):
    # plot mean distances given by function for different number of leaves and plot them
    mean_array = []
    var_array = []
    prgrs = 0.05
    for i in range(min_num_leaves,max_num_leaves):
        if (i/max_num_leaves >= prgrs):
            print('progress: ',prgrs)
            prgrs += 0.05
        if func == given_focal_tree_dist: # For now we assume that here we want to compute the distance from one caterpillar tree to num_tree_pairs coalescent trees
            statistics = func(i,num_tree_pairs, mean = True, focal_tree = sim.identity_caterpillar(i))
        else:
            statistics = func(i,num_tree_pairs, mean = True)
        diameter = (i-1)*(i-2)/2
        mean_array.append(statistics[0])
        var_array.append(statistics[1])
    if variance == True:
        d = pd.DataFrame(data=var_array)
        plts.plot_dots(d, [0,1], output_file, line = True)
    else:
        d = pd.DataFrame(data=mean_array, index = [i for i in range(min_num_leaves, max_num_leaves)])
        sns.scatterplot(data = d, Color = '#b02538', legend = False)
        plt.xlabel("Number of leaves")
        plt.ylabel("Mean distance")
        plt.savefig("../simulations/distance_distribution/coalescent/mean_dist_n_" + str(min_num_leaves) + "_to_" + str(max_num_leaves) + "_N_" + str(num_tree_pairs) + "_scatter.eps")
        # plt.show()
        plt.clf()
        sns.lineplot(data = d, Color = '#b02538', legend = False)
        plt.xlabel("Number of leaves")
        plt.ylabel("Mean distance")
        plt.savefig("../simulations/distance_distribution/coalescent/mean_dist_n_" + str(min_num_leaves) + "_to_" + str(max_num_leaves) + "_N_" + str(num_tree_pairs) + "_line.eps")
        # plt.show()
    # mean_diff = []
    # for i in range(0, len(mean_array)-1):
    #     mean_diff.append(mean_array[i+1] - mean_array[i])
    # d = pd.DataFrame(data=mean_diff)
    # plts.plot_dots(d, line = True)


def mean_distance_exp_n(func, max_exp, num_tree_pairs, output_file = ''):
    # plot mean distances given by function for different number of leaves and plot them
    mean_array = []
    # var_array = []
    for i in range(2,2 + max_exp):
        num_leaves = 2**i
        print(num_leaves, 'leaves')
        if func == given_focal_tree_dist: # For now we assume that here we want to compute the distance from one caterpillar tree to num_tree_pairs coalescent trees
            statistics = func(num_leaves,num_tree_pairs, mean = True, focal_tree = sim.identity_caterpillar(num_leaves))
        else:
            statistics = func(num_leaves,num_tree_pairs, mean = True)
        diameter = (num_leaves-1)*(num_leaves-2)/2
        mean_array.append(statistics[0])
        # var_array.append(statistics[1])
    print(mean_array)
    # print(var_array)
    d = pd.DataFrame(data=mean_array, index = [i for i in range(2, 2+max_exp)])
    sns.scatterplot(data = d, Color = '#b02538', legend = False)
    plt.xlabel("k for 2^k leaves")
    plt.ylabel("Mean distance")
    plt.savefig("../simulations/distance_distribution/coalescent/mean_dist_exp_n_k_2_to_" + str(max_exp) + "_N_" + str(num_tree_pairs) + "_scatter.eps")
    plt.show()
    # plt.plot(var_array)


def mean_distance_100_n(func, min_num_leaves, max_num_leaves, num_trees, output_file = ''):
    # plot mean distances given by function for different number of leaves and plot them
    mean_array = []
    # var_array = []
    for i in range(min_num_leaves,max_num_leaves):
        print(100*i)
        if func == given_focal_tree_dist: # For now we assume that here we want to compute the distance from one caterpillar tree to num_trees coalescent trees
            statistics = func(100*i,num_trees, mean = True, focal_tree = sim.identity_caterpillar(100*i))
        else:
            statistics = func(100 * i,num_trees, mean = True)
        diameter = (100 * i-1)*(100 *i-2)/2
        mean_array.append(statistics[0])
        # var_array.append(statistics[1])
    print(mean_array)
    # print(var_array)
    d = pd.DataFrame(data=mean_array)
    plts.plot_dots(d, [0,1], output_file)
    # plt.plot(var_array)


def mean_distance_repeat(func, num_leaves, num_iterations, num_trees, output_file = '', distance_file = '', balanced = False, variance = False):
    # plot mean distances given by function for different number of leaves and plot them
    diameter = (num_leaves-1)*(num_leaves-2)/2
    mean_array = []
    var_array = []
    for i in range(0,num_iterations):
        if balanced == True:
            statistics = given_focal_tree_dist(num_leaves, num_trees, focal_tree = sim.sim_cat(num_leaves, 1).trees[0], mean = True, distances_file=distance_file)
        else:
            statistics = func(num_leaves,num_trees, mean = True, distances_file=distance_file)
        mean_array.append(statistics[0])
        var_array.append(statistics[1])
    # print(mean_array)
    # print(var_array)
    if distance_file != '':
        np.savetxt(distance_file, mean_array)
    if variance == True:
        d = pd.DataFrame(data = var_array)
        print(np.var(var_array), min(var_array), max(var_array))
        plts.plot_dots(d)
    else:
        # d = pd.DataFrame(data =  [i/diameter for i in mean_array], columns = ["mean distance"]) # If distances in mean_array aren't normalised yet
        d = pd.DataFrame(data =  mean_array, columns = ["mean distance"])
        plts.plot_dots(d, [0,1], filehandle = output_file)


def expected_dist(num_leaves):
    # Returns a list of expected distances for trees on 3 to num_leaves trees
    if num_leaves == 2:
        return(0)
    exp_dist = [0]
    for n in range(3, num_leaves+1):
        p = []
        if n == 3:
            exp_dist[n-3] = 0
        else:
            exp_dist.append(exp_dist[n-4])
        for k in range(1, n):
            p.append(1) # values p_{k,n} for n=i
            for i in range(1, k):
                p[k-1] = p[k-1] * (1-2/((n-i+1)*(n-i)))
            p[k-1] = p[k-1] * (2/((n-k+1)*(n-k)))
            exp_dist[n-3] += p[k-1] * (k-1)
    return(exp_dist)


def plot_approx_exp_dist(max_num_leaves, output_file = ''):
    exp_dist = expected_dist(max_num_leaves)
    norm_exp_dist = []
    for i in range(3,max_num_leaves):
        diameter = (i-1)*(i-2)/2
        norm_exp_dist.append(exp_dist[i-3]/diameter)
    d = pd.DataFrame(data = norm_exp_dist)
    if output_file != '':
        plt.savefig(output_file)
    plts.plot_dots(d)


def lower_bound_expected_distance(num_leaves):
    # Compute a lower bound for the expected distance, using the average height of the tree and taking the mean distance of caterpillar trees with mx number of leaves that can be embedded according to average height -- This is a terrible lower bound, don't use it!
    # Compute harmonic number - 1:
    h = 0
    for i in range(2,num_leaves+1):
        h += 1/i
    return(2/3*(2*h**2 - h))


def compare_expected_dist_to_simulation(num_leaves, num_trees, output_file = '', all_elements = True):
    # Compare expected distances of trees on 3 to num_leaves leaves to mean of simulated distance distribution (based on num_trees simulated coalescent trees)
    # if all_elements == False: only do this for trees on num_leaves leaves.
    exp_dist = expected_dist(num_leaves)
    norm_exp_dist = []
    norm_sim_dist = []
    if all_elements == True:
        for i in range(3,num_leaves+1):
            diameter = (i-1)*(i-2)/2
            # print(Fraction(np.mean(exp_dist[i-3])).limit_denominator())
            norm_sim_dist.append(coal_pw_dist(i, num_trees, mean=True)[0])
            norm_exp_dist.append(exp_dist[i-3]/diameter)
    if all_elements == False:
        diameter = (num_leaves - 1)*(num_leaves - 2)/2
        norm_sim_dist.append(coal_pw_dist(num_leaves, num_trees, mean=True)[0])
        norm_exp_dist.append(exp_dist[num_leaves - 3]/diameter)
    d = pd.DataFrame(data = list(zip(norm_exp_dist, norm_sim_dist)), columns = ["approximated expectation", "mean of simulation"])
    sns.scatterplot(data=d, s = 50, legend = True)
    if output_file != '':
        plt.savefig(output_file)
    plt.show()


def plot_moves_per_iteration(num_leaves, num_trees, output_file = ''):
    diameter = (num_leaves - 1)*(num_leaves - 2)/2
    moves_per_iteration = []
    for i in range(0, num_leaves-1):
        moves_per_iteration.append(0)
    for k in range(0,num_trees):
        tree_list = sim.sim_coal(num_leaves, 2)
        fp_moves = findpath_moves_per_iteration(tree_list.trees[0], tree_list.trees[1])
        moves_this_iteration = []
        for i in range(0,num_leaves-1):
            moves_this_iteration.append(fp_moves[i])
        if sum(moves_this_iteration) != 0:
            moves_this_iteration = [i/sum(moves_this_iteration) for i in moves_this_iteration]
        for i in range(0,num_leaves-1):
            moves_per_iteration[i] += moves_this_iteration[i]
    for i in range(0,num_leaves-1):
        moves_per_iteration[i] = moves_per_iteration[i]
    # print(sum(moves_per_iteration))
    moves_per_iteration = [i/num_trees for i in moves_per_iteration]
    print(moves_per_iteration)
    print(sum(moves_this_iteration))
    d = pd.DataFrame(data = moves_per_iteration)
    # bins = np.arange(-.5, num_leaves + 0.5, 1)
    # plts.plot_hist(d, bins, density = False)
    plts.plot_dots(d, filehandle = output_file)
    # p = sns.scatterplot(data=d, s = 50)
    # plt.show()


def random_walk_distance(num_leaves, k, num_iterations, output_file = '', mean = False):
    # k: length of random walk
    distances = []
    diameter = (num_leaves - 1)*(num_leaves -2)/2
    for i in range(0, num_iterations):
        tree = sim.sim_coal(num_leaves,1).trees[0]
        distances.append(random_walk(tree, k))
    if mean == True:
        return(np.mean(distances))
    else:
        print(np.mean(distances))
        d = pd.DataFrame(data = distances)
        # bins = np.arange(-.5, k + 1.5, 1)
        # plts.plot_hist(d, bins, filehandle = output_file)
        sns.histplot(d, Color = '#b02538', Edgecolor = 'black', alpha=1, binwidth=1, binrange = [-.5,diameter+1.5], stat = 'density', legend = False)
        plt.xlabel("Distance")
        # plt.ylabel("Frequency")
        plt.ylabel("")
        plt.savefig("../simulations/distance_distribution/coalescent/random_walk_dist_" + str(num_leaves) + "_n_" + str(k) + "_steps_" + str(num_iterations) + "_repetitions_.eps")
        plt.show()

def random_walk_mean_distance(num_leaves, k_min, k_max, num_iterations, output_file = '', median = False):
    mean_dist = []
    for k in range(k_min, k_max):
        mean_dist.append(random_walk_distance(num_leaves, k, num_iterations, mean = True))
    if median == True:
        index = 0
        m = np.median(mean_dist)
        while(mean_dist[index] < m and index < len(mean_dist)):
            index += 1
        mean_mean_dist = mean_dist[index:]
        return(np.mean(mean_mean_dist), index) # Return approximation of value the distance converges to + length index of random walk from when it is closed convergence (Median)
        print(index)
    d = pd.DataFrame(data = mean_dist)
    # plts.plot_dots(d, filehandle = output_file, line = True)
    sns.scatterplot(data = d, Color = '#b02538', legend = False)
    plt.xlabel("Length of random walk")
    plt.ylabel("Mean distance")
    plt.savefig("../simulations/distance_distribution/coalescent/random_walk_mean_dist" + str(num_leaves) + "_n_" + str(k_min) + "_to_" + str(k_max) +"_k_" + str(num_iterations) +  "_num_repeats_scatter.eps")
    plt.show()

def random_walk_mean_distance_exp(num_leaves, k, num_iterations, output_file = '', median = False):
    mean_dist = []
    diameter = (num_leaves-1)*(num_leaves-2)/2
    for i in range(1, k+1):
        print(i)
        mean_dist.append(random_walk_distance(num_leaves, 2**i, num_iterations, mean = True)/diameter)
    if median == True:
        index = 0
        m = np.median(mean_dist)
        while(mean_dist[index] < m and index < len(mean_dist)):
            index += 1
        mean_mean_dist = mean_dist[index:]
        return(np.mean(mean_mean_dist), index) # Return approximation of value the distance converges to + length index of random walk from when it is closed convergence (Median)
        print(index)
    print(mean_dist)
    d = pd.DataFrame(data = mean_dist, index = [i for i in range(1, k+1)])
    # plts.plot_dots(d, filehandle = output_file, line = True)
    sns.scatterplot(data = d, Color = '#b02538', legend = False)
    plt.xlabel("k where 2^k is length of random walk")
    plt.ylabel("Mean distance")
    plt.savefig("../simulations/distance_distribution/coalescent/random_walk_exp_mean_dist" + str(num_leaves) + "_n_2_to_" + str(k) +"_k_" + str(num_iterations) +  "_num_repeats_scatter.eps")
    plt.show()


def number_rnni_edges(num_leaves):
    # Compute the number of RNNI edges in the graph on num_leaves leaves
    factorial = 1
    frac_sum = 0
    for i in range(1, num_leaves): #Compute (num_leaves-1)!
        factorial = factorial * i
        if i > 1:
            frac_sum = frac_sum + 1/i
    num_trees = num_leaves*factorial**2/2**(num_leaves-1)
    nni_edges = 2/3*num_trees*frac_sum
    rank_edges = 1/2 * num_trees * (num_leaves - 2 - 2*(frac_sum))
    return(3*nni_edges + rank_edges)


def expected_one_neighbourhood_size_n(num_leaves):
    neighbourhood_size = []
    for i in range(3, num_leaves):
        factorial = 1
        for j in range(1, i): #Compute (i-1)!
            factorial = factorial * j
        num_trees = i*factorial**2/2**(i-1)
        neighbourhood_size.append(1/num_trees * 2 * number_rnni_edges(i))
    d = pd.DataFrame(neighbourhood_size)
    plts.plot_dots(d, line = True)

if __name__ == '__main__':
    # pw_tree_list_dist("../../../../online/online-ms/non_bayes/output/trees/rooted_20_leaves_5_drops_100_repeats.trees", output_file="../../../../online/online-ms/non_bayes/output/trees/RNNI_rooted_20_leaves_5_drops_100_repeats.pdf", metric = 'RNNI', nexus = False)
    # coal_pw_dist_space_efficient(20,100000) # Output file!
    # print(number_rnni_edges(7))
    # for i in range(4, 40):
    #     diameter = (i-1)*(i-2)/2
    #     print(lower_bound_expected_distance(i)/diameter)
    # mean_distance_repeat(coal_pw_dist, 3, 40, 1000, variance=True)
    # compare_expected_dist_to_simulation(10000, 100, all_elements=False)
    # mean_distance_n(coal_pw_dist, 3, 10, 1000)
    # for num_leaves in range(3,10):
    #     print(random_walk_mean_distance(num_leaves,1,800,1000, median = True))#, output_file = '../simulations/distance_distribution/coalescent/random_walk_mean_dist_n_6_k_1_to_1000_N_1000.eps')
    random_walk_distance(10, 1000, 100000)#, output_file = '../simulations/distance_distribution/coalescent/random_walk_dist_n_6_k_20_N_1000.eps')
    # random_walk_mean_distance(10,1,4*36,1000)#, output_file = '../simulations/distance_distribution/coalescent/random_walk_mean_dist_n_6_k_1_to_1000_N_1000.eps')
    # random_walk_mean_distance_exp(10,15,10000)
    # random_walk_mean_distance(7,1,1000,1000, output_file = '../simulations/distance_distribution/coalescent/random_walk_mean_dist_n_7_k_1_to_1000_N_1000.eps')
    # coal_tree_list = sim.sim_coal(10, 10000)
    # tree_list = read_nexus('../simulations/posterior/coal/coal_alignment_20_sequences_10000_length.trees', ranked = True) # Count number of trees for posterior sample (simulated)
    # distances = rnni.pw_rnni_dist(tree_list)
    # count_equal_trees(tree_list) # Not saved anywhere, but this does (as expected) give a constant for uniform distribution. It could be used to investigate other distribution and show that we need a distance for further analyses

    # plot_approx_exp_dist(500, output_file = '../simulations/distance_distribution/coalescent/approx_exp_dist_n_3_to_500.eps')
    # plot_moves_per_iteration(100, 10000, output_file='../simulations/distance_distribution/coalescent/moves_per_iteration_n_100_N_10000.eps')
    # plot_moves_per_iteration(4, 100000)

    # compare_expected_dist_to_simulation(200, 1000)
    # compare_expected_dist_to_simulation(40, 20000, output_file='../simulations/distance_distribution/coalescent/compare_expected_dist_to_simulation_40_n_20000_N.eps')

    # compare_given_focal_tree_dist(16, 10000, focal_tree1 = sim.balanced_tree_16_leaves(), focal_tree2 = sim.sim_cat(16,1).trees[0]) #, output_file = '../simulations/distance_distribution/coalescent/compare_cat_balanced_16_n_10000_N.eps')
    # given_focal_tree_dist(16, 10000, sim.balanced_tree_16_leaves(), output_file = '../simulations/distance_distribution/coalescent/dist_distribution_to_fully_balanced_16_n_10000_N.eps')

    # dist_distribution_to_caterpillars(20,10000, output_file = '../simulations/distance_distribution/coalescent/dist_distribution_to_caterpillars_20_n_10000_N.eps')
    # dist_distribution_btw_caterpillars(20,20000, output_file = '../simulations/distance_distribution/coalescent/dist_distribution_btw_caterpillars_20_n_20000_N.eps')
    # mean_distance_n(dist_distribution_to_caterpillars, 3, 40, 10000, output_file = '../simulations/distance_distribution/coalescent/cat_mean_dist_n_3_to_40_N_10000.eps')
    # mean_distance_n(given_focal_tree_dist, 3, 40, 10000, output_file = '../simulations/distance_distribution/coalescent/from_cat_mean_dist_n_3_to_40_N_10000.eps')
    # mean_distance_n(dist_distribution_btw_caterpillars, 3, 40, 10000, output_file = '../simulations/distance_distribution/coalescent/btw_cat_mean_dist_n_3_to_40_N_10000.eps')
    # mean_distance_n(coal_pw_dist, 3, 300, 1000, output_file = '../simulations/distance_distribution/coalescent/mean_dist_n_3_to_300_N_1000.eps')
    # mean_distance_n(coal_pw_dist, 3, 100, 10000)
    # mean_distance_exp_n(coal_pw_dist_space_efficient, 10, 10000)
    # mean_distance_100_n(coal_pw_dist_space_efficient, 100, 101, 10) #, output_file = '../simulations/distance_distribution/coalescent/mean_dist_log_n_1_to_9_N_10000.eps')
    # [0.8285882293673883, 0.8300589130067364, 0.8301474363627227, 0.8307205053473854, 0.8311932083023579, 0.8310782928817221, 0.831468037175452, 0.8319265701059583, 0.8316929550362643, 0.8321324688541974, 0.8318900686265788, 0.8321025213961065, 0.8324846223507965, 0.8324120056764494, 0.8318717888630501, 0.8323381148416531]

    # mean_distance_repeat(dist_distribution_to_caterpillars, 16, 500, 10000, output_file = '../simulations/distance_distribution/coalescent/mean_caterpillar_distance_repeat_n_16_N_10000_500_iterations.eps', distance_file = '../simulations/distance_distribution/coalescent/mean_caterpillar_distance_repeat_n_16_N_10000_500_iterations.np')
    # mean_distance_repeat(given_focal_tree_dist, num_leaves = 16, num_iterations = 500, num_trees = 10000, output_file = '../simulations/distance_distribution/coalescent/mean_focal_distance_repeat_n_16_N_10000_500_iterations.eps', distance_file = '../simulations/distance_distribution/coalescent/mean_focal_distance_repeat_n_16_N_10000_500_iterations.np')
    # mean_distance_repeat(given_focal_tree_dist, 16, 500, 10000, output_file = '../simulations/distance_distribution/coalescent/mean_caterpillar_distance_repeat_n_16_N_10000_500_iterations.eps', distance_file = '../simulations/distance_distribution/coalescent/mean_from_caterpillar_distance_repeat_n_16_N_10000_500_iterations.np')
    # mean_distance_repeat(coal_pw_dist, 16, 500, 10000, output_file = '../simulations/distance_distribution/coalescent/mean_distance_repeat_n_16_N_10000_500_iterations.eps', distance_file = '../simulations/distance_distribution/coalescent/mean_distance_repeat_n_16_N_10000_500_iterations.np')
    # coal_pw_dist(20,20000, output_file = '../simulations/distance_distribution/coalescent/own_coal_distr_20_n_20000_N.eps')
    # caterpillar_dist_distribution(20,20000, output_file='../simulations/distance_distribution/coalescent/caterpillar_distances_20_n_20000_N.eps')
    # coal_focal_dist(20, 20000, output_file='../simulations/distance_distribution/coalescent/coal_focal_dist_20_n_20000_N.eps')

    # pw_tree_list_dist('../simulations/simulated_trees/coal/20000/coal_trees_20_n.nex', output_file = '../simulations/distance_distribution/coalescent/wrf_distribution_20_n_20000_N.eps', metric = 'wRF')
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

    # plts.mean_comparison()

    # mean_diff_list = []
    # for i in range(0, len(mean_dist_list)-1):
    #     mean_diff_list.append(mean_dist_list[i+1] - mean_dist_list[i])
    # plt.plot(mean_diff_list)
    # plt.show()
