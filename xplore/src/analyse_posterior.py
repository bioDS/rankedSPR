# Analyse posterior distributions
# Especially: Compare two distributions, using the RNNI distance
__author__ = 'Lena Collienne'

import os.path
import sys
# sys.path.append('../..')

from analyse_distance_distribution import *


def compare_two_samples(file1, file2, output_file = ''):
    # Compare the distance distributions of two samples separately with the distance distribution resulting from merging the two samples
    # Read tree_lists:
    tree_list1 = read_nexus(file1, ranked = True)
    tree_list2 = read_nexus(file2, ranked = True)
    # Merge tree lists to tree_list3:
    num_total_trees = tree_list1.num_trees + tree_list2.num_trees
    concat_trees = (TREE * num_total_trees)()
    for i in range(0, tree_list1.num_trees):
        concat_trees[i] = tree_list1.trees[i]
    for i in range(0, tree_list2.num_trees):
        concat_trees[tree_list1.num_trees+i] = tree_list2.trees[i]
    tree_list3 = TREE_LIST(concat_trees, num_total_trees)

    # Compute three distance matrices: tree_list1, 2, 3
    # For plotting we need num_trees and diameter:
    num_leaves = tree_list1.trees[0].num_leaves
    rnni_diameter = int((num_leaves-1)*(num_leaves-2)/2)

    # Plotting RNNI distances for all pairs T_i, T_j, i<j, output as list (saves a lot of memory + runtime for plot)
    print('Compute distances.')
    distances1 = rnni.pw_rnni_dist(tree_list1, list = True)
    distances2 = rnni.pw_rnni_dist(tree_list2, list = True)
    distances3 = rnni.pw_rnni_dist(tree_list3, list = True)
    print('Done computing distances.')

    # Plot histograms based on distances
    bins = np.arange(-.5, rnni_diameter + 1.5, 1)
    df = pd.DataFrame(data = list(zip(distances3, distances1, distances2)), columns = ["joint sample", "sample 1", "sample 2"])
    plts.plot_hist(df, bins, output_file)

def dist_to_mcc(summary_tree_file, tree_file):
    # Compute distance of all trees of sample to one tree (in summary_tree_file), both files nexus format!
    # summary_tree = read_nexus(summary_tree_file, ranked  = True).trees[0] #Careful: Leaves need to reveive same labels -- we might need to add summary tree to tree_file
    tree_list = read_nexus(tree_file, ranked = True)
    # summary_tree = tree_list.trees[tree_list.num_trees - 1]
    summary_tree = read_nexus(summary_tree_file, ranked = True).trees[0]
    distance = []
    # print(tree_list.num_trees)
    for i in range(0, tree_list.num_trees):
        distance.append(findpath_distance(summary_tree, tree_list.trees[i]))
    print(distance)
    plt.plot(distance)
    plt.show()

        



if __name__ == '__main__':
    dist_to_mcc("../simulations/posterior/primates/mcc.trees", "../simulations/posterior/primates/primates_small.trees")
    # trees = (TREE * 2)()
    # for i in range(1,3):
    #     with open('../simulations/posterior/comparison/tree_' + str(i) + '_on_20_leaves.new') as f:
    #         tree_string = f.read() # This assumes that the file actually only contains one line that is a tree in newick format
    #         trees[i-1] = read_newick(tree_string)
    # # for i in range(1,3):
    #     # all_pw_dist('../simulations/posterior/comparison/coal_alignment_tree_' + str(i) + '_on_20_sequences_10000_length.trees')
    #     # all_pw_dist('../simulations/posterior/comparison/coal_alignment_tree_1_on_20_sequences_10000_length.trees')
    # compare_two_samples('../simulations/posterior/comparison/coal_alignment_tree_1_on_20_sequences_10000_length.trees', '../simulations/posterior/comparison/coal_alignment_tree_2_on_20_sequences_10000_length.trees', output_file = '../simulations/posterior/comparison/coal_alignment_plot_distances.eps')