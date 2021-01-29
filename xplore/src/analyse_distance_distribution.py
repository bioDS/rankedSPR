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
import rf_distances as rf
import rnni_distances as rnni
import plots as plts


def all_pw_distances(input_file, output_file = '', matrix_file = '', metric = 'RNNI'):
    # plot and save histogram of all pw distances ('RNNI' od 'RF'). inpput_file is filehandle to file with input trees, matrix_file (if specified) is either file to read the matrix from (if already exists), or file to save matrix in (make sure it didn't exist before!), output_file is filehandle for saving the plot (histoogram).
    if metric == 'RNNI':
        # Read trees in C format (for RNNI distance computation)
        print("Read trees")
        tree_list = read_nexus(input_file, ranked = True)[0]
        print("Done reading trees")
        num_trees = tree_list.num_trees
        num_leaves = tree_list.trees[0].num_leaves
        rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

        # Plotting RF distances for all pairs T_i, T_j, i<j
        try:
            f = open(matrix_file)
            f.close()
            pw_distances_rnni = np.loadtxt(matrix_file, delimiter = ' ')
        except IOError:
            pw_distances_rnni = rnni.pw_rnni_dist(tree_list)
            if matrix_file != '':
                np.savetxt(matrix_file, pw_distances_rnni, delimiter = ' ')
        finally:
            plts.plot_hist(pw_distances_rnni, output_file)
            if output_file != '':
                plt.savefig(output_file)
    
    elif metric == 'RF':
        # Read trees in ete3 format:
        print("Read trees")
        tree_list, leaf_labels = rf.read_ete_nexus(input_file)
        print("Done reading trees")
        num_leaves = len(leaf_labels)
        rf_diameter = int(2*(num_leaves - 1))

        # Plotting RF distances for all pairs T_i, T_j, i<j
        try:
            f = open(matrix_file)
            f.close()
            pw_distances_rf = np.loadtxt(matrix_file, delimiter = ' ')
        except IOError:
            pw_distances_rf = rf.pw_rf_dist(tree_list, list = True)
            if matrix_file != '':
                np.savetxt(matrix_file, pw_distances_rf, delimiter = ' ')
        finally:
            plts.plot_hist(pw_distances_rf, output_file)
            if output_file != '':
                plt.savefig(output_file)


if __name__ == '__main__':
    
    all_pw_distances('../simulations/simulated_trees/coal/coal_trees_20_n_100_N.nex', '../simulations/simulated_trees/coal/output.eps', metric = 'RF')

    # plt.plot(rnni_mean_dist_n(100, 20000), linestyle = 'None', marker = 'o', markersize = 6) # coalescent
    # plt.plot(rnni_mean_dist_n(40, 20000, model = 'bd'), linestyle = 'None', marker = 'o', markersize = 6) # birth-death
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
    # tree_list, leaf_labels = rf.read_ete_nexus(filehandle)
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