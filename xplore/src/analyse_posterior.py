# Analyse posterior distributions
# Especially: Compare two distributions, using the RNNI distance
__author__ = 'Lena Collienne'

import os.path
import sys
import time

from analyse_distance_distribution import *
from rnni_distances import pw_rnni_dist


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
    plt.clf()
    bins = np.arange(-.5, rnni_diameter + 1.5, 1)
    df = pd.DataFrame(data = list(zip(distances3, distances1, distances2)), columns = ["joint sample", "sample 1", "sample 2"])
    plts.plot_hist(df, bins, output_file)

def dist_to_mcc(summary_tree_file, tree_file, output_file):

    start_time = time.perf_counter() # For displaying time of computation
    # Compute distance of all trees of sample to one tree (in summary_tree_file), both files nexus format!
    # summary_tree = read_nexus(summary_tree_file, ranked  = True).trees[0] #Careful: Leaves need to reveive same labels -- we might need to add summary tree to tree_file
    
    # tree_list = read_nexus(tree_file, ranked = True) # We cannot read in the whole tree file at once, as this needs too much space
    summary_tree = read_nexus(summary_tree_file, ranked = True).trees[0]
    # print(summary_tree.num_leaves)    
    distance = []

    # Read trees from tree_file line by line, as in read_nexus, but without storing all trees in one list.
    # Precompiled Regex for a line containing a tree
    re_tree = re.compile("\t?tree .*=? (.*);$", flags=re.I | re.MULTILINE)
    # Used to delete the ; and a potential branchlength of the root

    # running variables for reading trees
    i = 0
    # Regex to delete additional data in []
    # brackets = re.compile(r'\[[^\]]*\]')

    with open(tree_file, 'r') as f:
        # Read trees
        for line in f:
            if re_tree.match(line):
                # First extract the newick string and then delete everything after the last occurence of ')'
                tree_string = f'{re.split(re_tree, line)[1][:re.split(re_tree, line)[1].rfind(")")+1]};'
                # Delete data in [] from newick, otherwise read_newick breaks
                tree = read_newick(tree_string, ranked=True)
                distance.append(findpath_distance(summary_tree, tree))
                i+=1
                if (i%10000 == 0):
                    current_time = time.perf_counter()
                    print("tree number: ", i, f" time: {current_time - start_time:0.4f} seconds")
                    if (i==100000):
                        break

    end_time = time.perf_counter()
    print(f"time needed: {end_time - start_time:0.4f} seconds")
    plt.clf()
    plt.plot(distance)
    plt.savefig(output_file)


def visual_all_distances(tree_file, output_file):
    tree_list = read_nexus(tree_file, ranked = True)
    distances = pw_rnni_dist(tree_list, list=False)

    dist_list = pw_rnni_dist(tree_list, list=True)
    d = pd.DataFrame(data = dist_list)
    num_leaves = tree_list.trees[0].num_leaves
    rnni_diameter = (num_leaves-1)*(num_leaves-2)/2
    # sns.histplot(d, Color = '#b02538', Edgecolor = 'black', alpha=1, binwidth=1, binrange = [-.5,rnni_diameter+1.5], stat = 'density', legend = False)
    # plt.show()
    plt.clf()
    sns.heatmap(distances,vmax=5)
    plt.savefig(output_file)
    plt.show()


if __name__ == '__main__':
    # dist_to_mcc("../simulations/posterior/primates/mcc_6000.trees", "../simulations/posterior/primates/primates_normal.trees", "../simulations/posterior/primates/mcc_6000_dist.eps")
    consec_trees_dist("../simulations/posterior/primates/primates_normal.trees", "../simulations/posterior/primates/pw_dist_primates_normal.eps")
    # visual_all_distances("../simulations/posterior/primates/primates_small.trees", "../simulations/posterior/primates/all_dist_matrix_primates_small.eps")
    # trees = (TREE * 2)()
    # for i in range(1,3):
    #     with open('../simulations/posterior/comparison/tree_' + str(i) + '_on_20_leaves.new') as f:
    #         tree_string = f.read() # This assumes that the file actually only contains one line that is a tree in newick format
    #         trees[i-1] = read_newick(tree_string)
    # # for i in range(1,3):
    #     # all_pw_dist('../simulations/posterior/comparison/coal_alignment_tree_' + str(i) + '_on_20_sequences_10000_length.trees')
    #     # all_pw_dist('../simulations/posterior/comparison/coal_alignment_tree_1_on_20_sequences_10000_length.trees')
    # compare_two_samples('../simulations/posterior/comparison/coal_alignment_tree_1_on_20_sequences_10000_length.trees', '../simulations/posterior/comparison/coal_alignment_tree_2_on_20_sequences_10000_length.trees', output_file = '../simulations/posterior/comparison/coal_alignment_plot_distances.eps')