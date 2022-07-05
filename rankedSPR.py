__author__ = 'Lena Collienne'
# Computing the rankedSPR graph to test algorithms for computing distances for trees on a small number of leaves
import sys
sys.path.append('treeOclock/')
sys.path.append('treeOclock/dct_parser/')

import ctypes
import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import random
import re
from treeOclock.dct_parser.tree_io import *
from treeOclock import *
from simulate_trees import *
from os.path import exists
from numpy.ctypeslib import ndpointer


_seidel = ctypes.CDLL("./libseidel.so")
_seidel.test_function.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32)
_seidel.seidel.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32)
_seidel.seidel_recursive.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32, ctypes.c_int32)


# Compute the adjacency matrix of the rankedSPR graph
# If hspr = 1, compute matrix for rankedSPR graph, otherwise for HSPR graph
def rankedSPR_adjacency(num_leaves, hspr = 1):
    num_trees = math.factorial(num_leaves - 1) * math.factorial(num_leaves) / (2**(num_leaves - 1))
    tree_index = dict() # dict containing trees as keys (as strings of cluster representation) and their index in adjacency matrix as values.
    index = 0 # Index of the last added tree in tree_index
    not_visited = [] # Save the trees that are already generated, but have not been added to adjacency matrix, in list
    visited = [] # Save the trees that have already been used, i.e. their 1-neighbourhoods have already been considered
    start_tree = identity_caterpillar(num_leaves)
    start_tree_str = tree_to_cluster_string(start_tree)
    tree_index[start_tree_str] = 0
    current_tree = start_tree
    visited.append(start_tree_str)
    # not_visited.append(start_tree_str)
    index += 1

    # adjacency matrix, initialised to only have 0s:
    adj = np.zeros((int(num_trees), int(num_trees)))

    # Fill tree_index with all trees on num_leaves leaves by going through lots of one-neighbourhoods, starting at start_tree:
    # Note that this is VERY inefficient! It is only kind of quick for up to 7 leaves!
    # We could do this a lot more efficient if we used the coalescent process.
    while len(visited) < num_trees:
        neighbourhood = all_spr_neighbourhood(current_tree, hspr)
        current_tree_str = tree_to_cluster_string(current_tree)

        for i in range(0,neighbourhood.num_trees):
            tree = neighbourhood.trees[i]
            tree_str = tree_to_cluster_string(tree)
            if tree_str not in tree_index:
                tree_index[tree_str] = index
                # Add 1 to adjacency matrix:
                index += 1
            adj[tree_index[current_tree_str],tree_index[tree_str]]=1
            adj[tree_index[tree_str], tree_index[current_tree_str]]=1

        # Randomly take the tree for the next iteration from the list of neighbours
        next_tree = sim_coal(num_leaves,1).trees[0]
        # next_tree = neighbourhood.trees[next_tree_index]
        next_tree_str = tree_to_cluster_string(next_tree)
        # We might get stuck somewhere in space where we can only escape by picking a new (random) starting tree
        while (next_tree_str in visited):
            next_tree = sim_coal(num_leaves,1).trees[0] # randomly choose the next tree for which we will compute the one-neighbourhood
            next_tree_str = tree_to_cluster_string(next_tree)
        visited.append(next_tree_str)
        if next_tree_str not in tree_index:
            tree_index[next_tree_str] = index
            index += 1
        current_tree = next_tree # update current_tree

    # Save tree dict in file:
    # open file for writing
    if hspr == 0:
        f = open("SPR/tree_dict_" + str(num_leaves) + "_leaves_hspr.txt","w")
    else:
        f = open("SPR/tree_dict_" + str(num_leaves) + "_leaves.txt","w")

    # write file
    for key in tree_index:
        f.write(str(tree_index[key]) + " " +str(key))
        f.write("\n")

    # close file
    f.close()

    # Save adjacency matrix in file
    if hspr ==1:
        if not exists('SPR/adj_matrix_%s_leaves.npy' %num_leaves):
            np.save("SPR/adj_matrix_" + str(num_leaves) + "_leaves.npy", adj)
    else:
        if not exists('SPR/adj_matrix_%s_leaves_hspr.npy' %num_leaves):
            np.save("SPR/adj_matrix_" + str(num_leaves) + "_leaves_hspr.npy", adj)
    return(adj, tree_index)


def test_restricted_neighbourhood_search(num_leaves, num_tree_pairs, hspr = 1):
    # Compute adjacency matrix & distance matrix
    rspr_adj = rankedSPR_adjacency(num_leaves, hspr)
    rspr_distances = np.ascontiguousarray(rspr_adj[0], dtype=np.int32)
    _seidel.seidel(rspr_distances, rspr_distances.shape[0])

    # Now simulate trees to be used to check distance computation
    t_list = sim_coal(num_leaves, 2*num_tree_pairs)
    correct_distance = 0
    for i in range(0,num_tree_pairs):
        tree1_index = rspr_adj[1][tree_to_cluster_string(t_list.trees[i])]
        tree2_index = rspr_adj[1][tree_to_cluster_string(t_list.trees[i+1])]
        # print(tree_to_cluster_string(t_list.trees[i]))
        # print(tree_to_cluster_string(t_list.trees[i+1]))
        # print(rspr_distances[tree1_index][tree2_index], rankedspr_path_restricting_neighbourhood(t_list.trees[i],t_list.trees[i+1]))
        if (rspr_distances[tree1_index][tree2_index] == rankedspr_path_restricting_neighbourhood(t_list.trees[i],t_list.trees[i+1], hspr)):
            correct_distance += 1
        else:
            print("tree1:", tree_to_cluster_string(t_list.trees[i]))
            print("tree2:", tree_to_cluster_string(t_list.trees[i+1]))
            print("correct distance:", rspr_distances[tree1_index][tree2_index], "approximated distance:", rankedspr_path_restricting_neighbourhood(t_list.trees[i],t_list.trees[i+1],hspr))
    print('correct distance:', correct_distance, 'out of', num_tree_pairs)


def test_restricted_neighbourhood_search_caterpillar(num_leaves, num_trees, hspr = 1):
    # Compute adjacency matrix & distance matrix
    rspr_adj = rankedSPR_adjacency(num_leaves, hspr)
    rspr_distances = np.ascontiguousarray(rspr_adj[0], dtype=np.int32)
    _seidel.seidel(rspr_distances, rspr_distances.shape[0])

    # Now simulate trees to be used to check distance computation
    t_list = sim_cat(num_leaves, num_trees)
    ctree = identity_caterpillar(num_leaves)
    correct_distance = 0
    ctree_index = rspr_adj[1][tree_to_cluster_string(ctree)]
    for i in range(0,num_trees):
        tree_index = rspr_adj[1][tree_to_cluster_string(t_list.trees[i])]
        # print(tree_to_cluster_string(t_list.trees[i]))
        # print(tree_to_cluster_string(t_list.trees[i+1]))
        # print(rspr_distances[tree_index][tree2_index], rankedspr_path_restricting_neighbourhood(t_list.trees[i],t_list.trees[i+1]))
        if (rspr_distances[tree_index][ctree_index] == rankedspr_path_restricting_neighbourhood(t_list.trees[i], ctree, hspr)):
            correct_distance += 1
        else:
            print("tree1:", tree_to_cluster_string(t_list.trees[i]))
            print("tree2:", tree_to_cluster_string(ctree))
            print("correct distance:", rspr_distances[tree_index][ctree_index], "approximated distance:", rankedspr_path_restricting_neighbourhood(t_list.trees[i],ctree,hspr))
    print('correct distance:', correct_distance, 'out of', num_trees)


def test_top_down_neighbourhood_search(num_leaves, num_tree_pairs):
    # Compute adjacency matrix & distance matrix
    rspr_adj = rankedSPR_adjacency(num_leaves)
    rspr_distances = np.ascontiguousarray(rspr_adj[0], dtype=np.int32)
    _seidel.seidel(rspr_distances, rspr_distances.shape[0])

    # Now simulate trees to be used to check distance computation
    t_list = sim_coal(num_leaves, 2*num_tree_pairs)
    correct_distance = 0
    for i in range(0,num_tree_pairs):
        tree1_index = rspr_adj[1][tree_to_cluster_string(t_list.trees[i])]
        tree2_index = rspr_adj[1][tree_to_cluster_string(t_list.trees[i+1])]
        # print("tree pair:")
        # print(tree_to_cluster_string(t_list.trees[i]))
        # print(tree_to_cluster_string(t_list.trees[i+1]))
        # print(rspr_distances[tree1_index][tree2_index], rankedspr_path_restricting_neighbourhood(t_list.trees[i],t_list.trees[i+1]))
        if (rspr_distances[tree1_index][tree2_index] == rankedspr_path_top_down_symm_diff(t_list.trees[i],t_list.trees[i+1])):
            correct_distance += 1
        else:
            print("tree1:", tree_to_cluster_string(t_list.trees[i]))
            print("tree2:", tree_to_cluster_string(t_list.trees[i+1]))
            print("correct distance:", rspr_distances[tree1_index][tree2_index], "approximated distance:", rankedspr_path_top_down_symm_diff(t_list.trees[i],t_list.trees[i+1]))
    print('correct distance:', correct_distance, 'out of', num_tree_pairs)


# Very slow and inefficient implementation of BFS for rankedSPR -- only useful for VERY small number of leaves
def rankedspr_bfs(start_tree, dest_tree, hspr=1):
    num_leaves = start_tree.num_leaves
    tree_dict = dict() # save trees (as cluster strings) and an index for each tree as value, so we can recover the path after running BFS (backtracking)
    index_dict = dict() # reverse of tree_dict (indices as keys and trees as values)
    predecessor = []
    to_visit = [] # queue containing next trees to be visited in BFS

    dest_tree_string = tree_to_cluster_string(dest_tree)
    
    # Initialise path?
    current_tree = start_tree

    tree_dict[tree_to_cluster_string(start_tree)] = 0
    index_dict[0] = tree_to_cluster_string(start_tree)
    index = 1 # index of the tree we currently consider (to be put as value for that tree into tree_dict)
    to_visit.append(current_tree)
    found = False # True if we found dest_tree
    # Start BFS
    while found == False:
        current_tree = to_visit.pop(0)
        current_tree_str = tree_to_cluster_string(current_tree)
        neighbours = all_spr_neighbourhood(current_tree,hspr)
        for i in range(0,neighbours.num_trees):
            tree = neighbours.trees[i]
            neighbour_string = tree_to_cluster_string(tree)
            if neighbour_string not in tree_dict:
                to_visit.append(tree)
                tree_dict[neighbour_string] = index
                index_dict[index]=neighbour_string
                predecessor.append(tree_dict[current_tree_str])
                index+=1
            if neighbour_string == dest_tree_string:
                found = True
                break

    # backtracking
    current_index = tree_dict[tree_to_cluster_string(dest_tree)]
    path_indices = [current_index]
    while (predecessor[current_index-1] != 0):
        path_indices.append(predecessor[current_index-1])
        current_index = predecessor[current_index-1]
    path_indices.append(0)
    # now turn path_indices array into path:
    path = []
    for i in range(len(path_indices)-1, -1, -1):
        path.append(index_dict[path_indices[i]])
    return(path)


# Find the ranks on which moves are performed on shortest path resulting from BFS
def bfs_path_rank_sequence(tree1, tree2):
    path = rankedspr_bfs(tree1, tree2, hspr=0)
    rank_list = []
    rank_count = []
    for i in range(0,len(path)-1):
        # for each move, find the lowest rank for which induced cluster changes -- this is the rank on which the HSPR move happened
        path[i] = str(path[i])
        path[i+1] = str(path[i+1])
        rank = 0
        for j in range(0, len(path[i])-1):
            if path[i][j] == '{':
                rank += 1
            if path[i][j] != path[i+1][j]:
                rank_list.append(rank)
                break
    for i in range(1,tree1.num_leaves-1):
        rank_count.append(rank_list.count(i))
    if max(rank_count) > 3:
        print("There is a rank for which more than one move is needed. The corresponding path is:")
        for tree in path:
            print(tree)
    return(rank_count)


def check_HSPR_moves_per_rank(num_leaves, num_tree_pairs):
    # simulate num_tree_pairs trees and check how moves are distributed across ranks in the trees on shortest path computed by bfs
    for i in range(0,num_tree_pairs):
        tree_list = sim_coal(num_leaves,2) # Simulate a pair of trees instead
        rank_count = bfs_path_rank_sequence(tree_list.trees[0], tree_list.trees[1])
        print(rank_count)


# use own implementation of coalescent to plot ranked SPR distances (HSPR if hspr=0, otherwise (default) RSPR) between coalescent trees (i.e. uniform ranked trees)
def coal_pw_spr_dist(num_leaves, num_tree_pairs, hspr = 1, output_file = '', distances_file = ''):
    # Plotting the distances for num_tree_pairs simulated pairs of trees and save plot (if filehandle given) in output_file
    distances = []

    for i in range(0,int(num_tree_pairs)):
        if i%100 == 0:
            print('iteration', i)
        tree_list = sim_coal(num_leaves,2) # Simulate a pair of trees instead of a list with num_tree trees
        distances.append(len(rankedspr_bfs(tree_list.trees[0], tree_list.trees[1]))-1)
    if distances_file != '':
        np.savetxt(distances_file,  distances, delimiter = ' ')
    # Plot histogram
    d = pd.DataFrame(data=distances)
    upper_bound = max(distances)
    b = np.arange(-.5, upper_bound + 1.5, 1)
    sns.set_theme(font_scale=1.2)
    sns.histplot(d, color = '#b02538', edgecolor = 'black', alpha=1, binwidth=1, binrange = [-.5,upper_bound+1.5], stat = 'density', legend = False)
    plt.xlabel("Distance")
    plt.ylabel("Proportion of trees")
    if hspr == 1:
        plt.savefig("SPR/plots/rspr_distribution_" + str(num_leaves) + "_n_" + str(num_tree_pairs) + ".eps")
    else:
        plt.savefig("SPR/plots/hspr_distribution_" + str(num_leaves) + "_n_" + str(num_tree_pairs) + ".eps")
    plt.clf()
    # plt.show()
    # plts.plot_hist(distances, bins, output_file)


# use own implementation of coalescent to plot difference in RSPR/HSPR distances between two coalescent trees and the same trees with one leaf deleted
def distance_del_leaf(num_leaves, num_deletions, num_tree_pairs, hspr = 1, output_file = '', distances_file = ''):
    # Plotting the distances for num_tree_pairs simulated pairs of trees and save plot (if filehandle given) in output_file
    distances = []

    for i in range(0,int(num_tree_pairs)):
        if i%100 == 0:
            print('iteration', i)
        tree_list = sim_coal(num_leaves,2) # Simulate a pair of trees instead of a list with num_tree trees
        tree1 = tree_list.trees[0]
        tree2 = tree_list.trees[1]
        d = len(rankedspr_bfs(tree_list.trees[0], tree_list.trees[1]))-1
        # take trees that result from deleting the same (randomly chosen) leaf
        # for i in range(0,num_deletions):
        #     r = random.randint(0,num_leaves-1-num_deletions)
        #     tree1 = del_leaf(tree1,r)
        #     tree2 = del_leaf(tree2,r)
        # d1 = len(rankedspr_bfs(tree1, tree2))-1
        # distances.append(d-d1)

        # # alternatively: try to delete every pair of leaves and look at minimum distance
        # current_dist = []
        # for i in range(0,num_leaves-1):
        #     tree1 = del_leaf(tree_list.trees[0],i)
        #     tree2 = del_leaf(tree_list.trees[1],i)
        #     for j in range(0,num_leaves-2):
        #         tree1_1 = del_leaf(tree1,j)
        #         tree2_1 = del_leaf(tree2,j)
        #         current_dist.append(len(rankedspr_bfs(tree1_1, tree2_1))-1)
        # distances.append(d - max(current_dist))

        # even another alternative: delete the two cherry leaves
        c1 = min(tree1.tree[num_leaves].children[0], tree1.tree[num_leaves].children[1])
        c2 = max(tree1.tree[num_leaves].children[0], tree1.tree[num_leaves].children[1])
        tree1 = del_leaf(tree1, c2)
        tree1 = del_leaf(tree1, c1)
        tree2 = del_leaf(tree2, c2)
        tree2 = del_leaf(tree2, c1)
        d1 = len(rankedspr_bfs(tree1, tree2))-1
        distances.append(d-d1)

        # if d-d1 == 3:
        #     print("original trees:")
        #     print(tree_to_cluster_string(tree_list.trees[0]))
        #     print(tree_to_cluster_string(tree_list.trees[1]))
        #     print("trees after deleting leaves:")
        #     print(tree_to_cluster_string(tree1))
        #     print(tree_to_cluster_string(tree2))

    print("maximum distance:", max(distances))
    if distances_file != '':
        np.savetxt(distances_file,  distances, delimiter = ' ')
    # Plot histogram
    d = pd.DataFrame(data=distances)
    upper_bound = max(distances)
    b = np.arange(-.5, upper_bound + 1.5, 1)
    sns.set_theme(font_scale=1.2)
    sns.histplot(d, color = '#b02538', edgecolor = 'black', alpha=1, binwidth=1, binrange = [-.5,upper_bound+1.5], stat = 'density', legend = False)
    plt.xlabel("Distance")
    plt.ylabel("Proportion of trees")
    if hspr == 1:
        plt.savefig("SPR/plots/rspr_dist_diff_" + str(num_leaves) + "_n_" + str(num_tree_pairs) + ".eps")
    else:
        plt.savefig("SPR/plots/hspr_dist_diff_" + str(num_leaves) + "_n_" + str(num_tree_pairs) + ".eps")
    plt.clf()
    # plt.show()
    # plts.plot_hist(distances, bins, output_file)


# use own implementation of coalescent to compare RSPR and HSPR distances between trees drawn from uniform distribution
def compare_hspr_rspr_uniform(num_leaves, num_tree_pairs, distances_file = ''):
    # Plotting the distances for num_tree_pairs simulated pairs of trees and save distances (if file handle provided)
    distances = [] # contains HSPR-RSPR distance for all simulated tree pairs
    for i in range(0,int(num_tree_pairs)):
        if i%100 == 0:
            print('iteration', i)
        tree_list = sim_coal(num_leaves,2) # Simulate a pair of trees instead of a list with num_tree trees
        distances.append(len(rankedspr_bfs(tree_list.trees[0], tree_list.trees[1], hspr=0))-len(rankedspr_bfs(tree_list.trees[0], tree_list.trees[1],hspr=1)))
    if distances_file != '':
        np.savetxt(distances_file,  distances, delimiter = ' ')
    # Plot histogram
    d = pd.DataFrame(data=distances)
    upper_bound = max(distances)
    b = np.arange(-.5, upper_bound + 1.5, 1)
    sns.set_theme(font_scale=1.2)
    sns.histplot(d, color = '#b02538', edgecolor = 'black', alpha=1, binwidth=1, binrange = [-.5,upper_bound+1.5], stat = 'density', legend = False)
    plt.xlabel("Distance")
    plt.ylabel("Proportion of trees")
    plt.savefig("SPR/plots/rspr_hspr_difference_" + str(num_leaves) + "_n_" + str(num_tree_pairs) + ".eps")
    plt.show()
    plt.clf()


# use own implementation of coalescent to compare differences between RSPR and HSPR shortest paths -- might be useful to determine if rank moves are always at beginning or end of shortest paths
def compare_hspr_rspr(num_leaves, num_tree_pairs):
    for i in range(0,int(num_tree_pairs)):
        if i%100 == 0:
            print('iteration', i)
        tree_list = sim_coal(num_leaves,2) # Simulate a pair of trees instead of a list with num_tree trees
        hspr_path = rankedspr_bfs(tree_list.trees[0], tree_list.trees[1], hspr=0)
        rspr_path = rankedspr_bfs(tree_list.trees[0], tree_list.trees[1], hspr=1)
        if(len(hspr_path)!=len(rspr_path)):
            print("hspr distance:", len(hspr_path)-1, "rspr_dist:", len(rspr_path)-1)
            print("hspr path:")
            for tree in hspr_path:
                print(tree)
            print("rspr path:")
            for tree in rspr_path:
                print(tree)


def find_rank_moves(num_leaves, num_tree_pairs):
    # simulate num_tree_pairs pairs of trees on num_leaves leaves and check the position of rank moves on shortest paths between those trees in RSPR (using BFS)
    # prints shortest paths that have rank moves not at the start or beginning, but somewhere in the middle
    for i in range(0,num_tree_pairs):
        tree_pair = sim_coal(num_leaves, 2)
        path = rankedspr_bfs(tree_pair.trees[0], tree_pair.trees[1])
        for i in range(0,len(path)-1):
            tree1_sorted = ''.join(sorted(str(path[i])))
            tree2_sorted = ''.join(sorted(str(path[i+1])))
            if tree1_sorted == tree2_sorted and i != 0 and i != len(path)-2:
                # two trees are connected by a rank move if the cluster strings contain exactly the same characters (are permutations of each other)
                print("rank move between:")
                print(path[i])
                print(path[i+1])
                print("entire path:")
                for tree in path:
                    print(tree)


def caterpillar_diameter_trees(n, hspr=1):
    # Checking which trees have diameter distance from identity caterpillar
    print("Reading trees")
    if hspr ==1:
        file = open('SPR/tree_dict_' + str(n) + '_leaves.txt')
        d = np.load('SPR/distance_matrix_' + str(n) + '_leaves.npy')
    else:
        file = open('SPR/tree_dict_' + str(n) + '_leaves_hspr.txt')
        d = np.load('SPR/distance_matrix_' + str(n) + '_leaves_hspr.npy')
    content = file.readlines()
    print("Done reading trees")

    max_indices = np.where(d == np.amax(d))
    max_coordinates = list(zip(max_indices[0], max_indices[1]))
    print("number of tree pair with diameter distance:")
    print(len(max_coordinates)/2)
    count = 0
    num_max_dist = 0
    # print(len(content[0]), content[0])
    for index1 in max_indices[1]:
        if max_indices[0][count] != 0:
            break
        count += 1
        print(index1, len(content[index1]), content[index1])
        # We only need to compare against first tree in file, as this is identity caterpillar tree (bc of symmetry)
        # if content[index1].count(',') == content[0].count(','): # check if content[index2] is caterpillar tree
        #     # print(content[0], content[index1])
        #     num_max_dist +=1

    print("number of caterpillar trees with diameter distance from identity caterpillar:", num_max_dist)
    print("total number of trees with diameter distance from identity caterpillar:", count)


def orbit_sizes(n, hspr=1):
    # find the number of trees at distance k from any tree in the distance matrix computed by SEIDEL
    # Output is an array of orbit sizes, where unique ones are only given once (e.g. all orbit sizes for same ranked topology will be the same)

    num_trees = int(math.factorial(n) * math.factorial(n-1) / (2**(n-1)))
    print("Start reading distance matrix")
    if hspr == 1:
        d = np.load('SPR/distance_matrix_' + str(n) + '_leaves.npy')
    else:
        d = np.load('SPR/distance_matrix_' + str(n) + '_leaves_hspr.npy')
    print("Done reading distance matrix")
    orbit_size = np.zeros((int(num_trees), int(np.amax(d)+1))) # initialise orbit sizes as zero matrix

    for i in range(0,num_trees):
        for j in range(0,np.amax(d)+1):
            orbit_size[i][j] = np.count_nonzero(d[i]==j)
            # print("distance", i, ":", num_trees, "trees")
    unique_rows = np.unique(orbit_size, axis=0, return_index = True)
    return(unique_rows) # unique_rows[0] contains the unique orbit sizes and unique_rows[1] contains the indices belonging to trees having those orbit sizes.

def print_orbits_with_trees(n, hspr=1):
    orbits = orbit_sizes(n, hspr)
    if hspr == 0:
        f = open('SPR/tree_dict_' + str(n) + '_leaves_hspr.txt')
    else:
        f = open('SPR/tree_dict_' + str(n) + '_leaves.txt')
    trees = f.readlines()
    for i in range(0,len(orbits[0])):
        print(str(trees[orbits[1][i]]), orbits[0][i])
    f.close()


def orbit_count_repetitions(tree, hspr=1):
    # print for every tree in 2-NH of tree how often it is counted if we perform two HSPR/RSPR moves on tree
    one_orbit_dict = dict()
    two_orbit_dict = dict()
    if hspr == 0:
        one_orbit = hspr_neighbourhood(tree)
    else:
        one_orbit = spr_neighbourhood(tree)
    # add all trees in 1-NH to one_orbit_dict
    for i in range(0,one_orbit.num_trees):
        one_orbit_dict[tree_to_cluster_string(one_orbit.trees[i])] = 1
    for i in range(0,one_orbit.num_trees):
        if hspr == 0:
            two_orbit = hspr_neighbourhood(one_orbit.trees[i])
        else:
            two_orbit = spr_neighbourhood(one_orbit.trees[i])
        for j in range(0,two_orbit.num_trees):
            tree_str = tree_to_cluster_string(two_orbit.trees[j])
            if tree_str in one_orbit_dict:
                one_orbit_dict[tree_str] += 1
            elif tree_str in two_orbit_dict:
                two_orbit_dict[tree_str] += 1
            elif tree_str == tree_to_cluster_string(tree): # ignore the initial tree
                pass
            else:
                two_orbit_dict[tree_str] = 1
    print("Size of 2-orbit:", len(two_orbit_dict))
    print("Size of 2-NH:", 1+len(one_orbit_dict)+len(two_orbit_dict))
    return(one_orbit_dict, two_orbit_dict)


# Check how many shortest tree pairs have shortest paths with a caterpillar tree on them (only possible for a small number of leaves!)
def check_caterpillar_on_shortest_path(num_leaves, hspr=1):
    if hspr == 1:
        d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves.npy')
    else:
        d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves_hspr.npy')

    if hspr == 0:
        f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves_hspr.txt', 'r')
    else:
        f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves.txt', 'r')

    num_ctrees_on_paths = 0

    num_trees = len(d)

    c_tree = identity_caterpillar(num_leaves)
    c_indices = [] # list of indices in tree_dict that correspond to caterpillar trees
    count = 0
    for line in f:
        tree_str = line[line.rfind('b')+1:]
        tree_str = tree_str[:tree_str.rfind(' ')]
        if same_topology(read_from_cluster(tree_str), c_tree)==0:
            c_indices.append(count)
        count += 1
    f.close()

    for i in range(0,num_trees):
        if math.floor(i%(num_trees/100)) == 0:
            print("progress:", i/num_trees)
        for j in range(i+1, num_trees):
            flag = False
            for k in c_indices: # check for every tree k if it is on a shortest i-j-path
                if d[i][k] + d[k][j] == d[i][j]:
                    flag = True
                    num_ctrees_on_paths += 1
                    break
                if flag == True:
                    break
            # if flag==False:
            #     print(i, j)
    print(num_ctrees_on_paths, "out of", (num_trees**2-num_trees)/2, "tree pairs have a path that has at least one caterpillar tree in them")
    return(num_ctrees_on_paths, num_trees)


def print_trees_at_diameter(num_leaves, hspr=1):
    if hspr == 1:
        d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves.npy')
        f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves.txt', 'r')
    else:
        d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves_hspr.npy')
        f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves_hspr.txt', 'r')
    tree_dict = f.readlines()
    d_max = np.amax(d)
    # print('trees at diameter distance:')
    tree_pairs = [] # actual trees at diameter distance (only one per topology pair)
    count = 0
    for coord in np.argwhere(d == d_max):
        if count%100==0:
            print("number of tree pairs considered: ", count)
        count+=1
        tree1_str = tree_dict[coord[0]].split("'")[1]
        tree2_str = tree_dict[coord[1]].split("'")[1]
        tree1 = read_from_cluster(tree1_str)
        tree2 = read_from_cluster(tree2_str)
        # print('tree1:', tree1_str)
        # print('tree2:', tree2_str)
        if len(tree_pairs) == 0:
            tree_pairs.append(set([tree1_str, tree2_str]))
        # print('tree pairs:', tree_pairs)
        topology_already_counted = False
        for i in range(0,len(tree_pairs)):
            tp = tree_pairs[i]
            # print('tp:', tp)
            t1_str = tp.pop()
            t2_str = tp.pop()
            t1 = read_from_cluster(t1_str)
            t2 = read_from_cluster(t2_str)
            tree_pairs[i] = set([t1_str, t2_str])
            if ((same_topology(t1,tree1)+same_topology(t2,tree2)== 0) or (same_topology(t2,tree1)+same_topology(t1,tree2)==0)):
                topology_already_counted = True
                break
        if topology_already_counted == False:
            tree_pairs.append(set([tree1_str, tree2_str]))
    print('tree pairs at diameter distance:')
    for pair in tree_pairs:
        print(pair)
    # print("diameter dist trees per topology:")
    # for topology in all_topologies:
    #     count = 0
    #     for pair in topology_pairs:
    #         if topology in pair:
    #             count +=1
    #     print("number of neighbours for topology", topology, ":", count)

    print("number of tree topology pairs:", len(tree_pairs))


# check if two trees have same unranked topology]
def same_unranked_tree(tree1, tree2):
    cluster_pattern = r'\{[^\}]*\}'
    tree1_str = str(tree_to_cluster_string(tree1))
    clusters1 = re.findall(cluster_pattern, tree1_str)
    tree2_str = str(tree_to_cluster_string(tree2))
    clusters2 = re.findall(cluster_pattern, tree2_str)
    if sorted(clusters1) == sorted(clusters2):
        return(True)
    else:
        return(False)


# Find longest shortest paths with only rank moves on them
def longest_rank_shortest_path(num_leaves):
    d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves.npy')
    f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves.txt', 'r')

    max_dist = np.amax(d)
    print('Diameter:', max_dist)

    tree_strings = f.readlines()

    # Put all trees into a dict (note that indices are sorted increasingly in file)
    index = 0
    tree_dict = dict()
    tree_index_dict = dict()
    for tree_str in tree_strings:
        tree_str = tree_str.split("'")[1]
        tree_dict[tree_str]=index
        tree_index_dict[index]=tree_str
        index += 1

    current_d = max_dist
    found_path = False # did we find a path with only rank moves on it?
    while(found_path==False):
        print("current distance:", current_d)
        print("number of trees at this distance:", len(np.argwhere(d==current_d)))
        for coord in np.argwhere(d == current_d):
            tree1_str = tree_index_dict[coord[0]]
            tree2_str = tree_index_dict[coord[1]]
            tree1 = read_from_cluster(tree1_str)
            tree2 = read_from_cluster(tree2_str)
            if same_unranked_tree(tree1, tree2):
                # Get length of shortest path that only has rank moves
                num_rank_moves = shortest_rank_path(tree1,tree2)
                if num_rank_moves == current_d:
                    found_path = True
                    print('Maximum length of a rank move only path is:', current_d)
                    print('Given by trees:\n', tree1_str, "\n", tree2_str)
                    return(0)
        current_d -=1


# simulate two trees (coalescent) and see whether there is a tree in 1-NH of starting tree resulting from rank move
def path_rank_moves_first(num_leaves, num_repeats):
    d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves.npy')
    f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves.txt', 'r')

    # Put all trees into a dict (note that indices are sorted increasingly in file)
    tree_strings = f.readlines()
    index = 0
    tree_dict = dict()
    tree_index_dict = dict()
    for tree_str in tree_strings:
        tree_str = tree_str.split("'")[1]
        tree_dict[tree_str]=index
        tree_index_dict[index]=tree_str
        index += 1

    for j in range(0,num_repeats):
        tree_list = sim_coal(num_leaves,2)
        tree1 = tree_list.trees[0]
        tree2 = tree_list.trees[1]
        tree1_str = str(tree_to_cluster_string(tree1)).split("'")[1]
        tree2_str = str(tree_to_cluster_string(tree2)).split("'")[1]

        tree1_index = tree_dict[tree1_str]
        tree2_index = tree_dict[tree2_str]

        rank_neighbours = all_rank_neighbours(tree1)
        for i in range(0,rank_neighbours.num_trees):
            neighbour_str = str(tree_to_cluster_string(rank_neighbours.trees[i])).split("'")[1]
            neighbour_index = tree_dict[neighbour_str]
            if d[tree1_index][neighbour_index] + d[neighbour_index][tree2_index] == d[tree1_index][tree2_index]:
                print("start:", tree1_str)
                print("neighbour:", neighbour_str)
                print("end:", tree2_str)


# Compute the maximum number of rank moves on a shortest path in RSPR (using the distance matrix for the whole tree space computed by SEIDEL)
def max_rank_move_shortest_path(tree1, tree2):
    num_leaves = tree1.num_leaves
    d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves.npy')
    f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves.txt', 'r')

    # Put all trees into a dict (note that indices are sorted increasingly in file)
    tree_strings = f.readlines()
    index = 0
    tree_dict = dict()
    tree_index_dict = dict()
    for tree_str in tree_strings:
        tree_str = tree_str.split("'")[1]
        tree_dict[tree_str]=index
        tree_index_dict[index]=tree_str
        index += 1

    tree1_str = str(tree_to_cluster_string(tree1)).split("'")[1]
    tree2_str = str(tree_to_cluster_string(tree2)).split("'")[1]

    tree1_index = tree_dict[tree1_str]
    tree2_index = tree_dict[tree2_str]

    distance = d[tree1_index][tree2_index]

    # for every tree that is on a shortest path, save all predecessors of it in dictionary pred:
    pred = dict()
    for tree_index in range(0,len(d)):
        if d[tree1_index][tree_index] + d[tree_index][tree2_index] == distance:
            tree = read_from_cluster(tree_index_dict[tree_index])
            neighbourhood = spr_neighbourhood(tree)
            for i in range(0, neighbourhood.num_trees):
                predecessor = neighbourhood.trees[i]
                pred_str = str(tree_to_cluster_string(predecessor)).split("'")[1]
                pred_index = tree_dict[pred_str]
                if d[tree1_index][pred_index] + d[pred_index][tree_index] + d[tree_index][tree2_index] == distance: # if predecessor is on shortest path from tree1 to tree2
                    if tree_index in pred:
                        pred[tree_index].add(pred_index)
                    else:
                        pred[tree_index] = set([pred_index])

    # print("starting tree:", tree1_str)
    # print("destination tree:", tree2_str)
    # for i in pred:
    #     for j in pred[i]:
    #         print(tree_index_dict[j], tree_index_dict[i])
    # print(len(pred))

    max_rank_moves = 0

    # print("start second part")
    # print(pred)

    # We now need to transform the predecessor dict into actual shortest paths and count how many rank moves are on each of these paths.
    found = False
    while True:
        current_path_rank_moves = 0
        # build path from end to beginning. Delete trees from pred dict, if all pred (i.e. all shortest path containing that tree) are considered.
        last_tree_index = tree2_index
        last_tree = read_from_cluster(tree_index_dict[tree2_index])
        last_popped = -1 # index of the last tree removed from pred[last_pooped_pred]
        last_popped_pred = -1
        while last_tree_index != tree1_index:
            # print(pred[last_tree_index])
            tree_index = pred[last_tree_index].pop()
            # print(tree_index)

            # we only delete the index out of the pred values if it was the last one where there were multiple choices on the currently computed path.
            if len(pred[last_tree_index]) >=1:
                if last_popped != -1:
                    if last_popped_pred in pred:
                        pred[last_popped_pred].add(last_popped)
                    else:
                        pred[last_popped_pred] = set([last_popped])
                last_popped = tree_index
                last_popped_pred = last_tree_index
            else:
                pred[last_tree_index].add(tree_index)
            # if tree_index in pred and len(pred[tree_index])>1: # if there are further paths going through tree_index, we add it back to the predecessor list of last_tree_index
            #     pred[last_tree_index].add(tree_index)
            if len(pred[last_tree_index])==0: # delete empty sets from pred (all paths through corresponding tree have already been considered)
                pred.pop(last_tree_index)
            tree_str = tree_index_dict[tree_index]
            tree = read_from_cluster(tree_str)
            if same_unranked_tree(tree,last_tree):
                current_path_rank_moves += 1
            # else:
            #     print("HSPR move")
            last_tree_index = tree_index
            last_tree = tree
        if current_path_rank_moves > max_rank_moves:
            max_rank_moves = current_path_rank_moves
        # print(current_path_rank_moves)
        if found == True:
            break
        done = True
        for i in pred:
            if len(pred[i]) != 1:
                done = False
        if done == True:
            found = True
    return(max_rank_moves)


# check max number of rank moves on shortest path for all tree pairs in RSPR
def check_max_rank_move_shortest_path(num_leaves):
    f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves.txt', 'r')

    # Put all trees into a dict (note that indices are sorted increasingly in file)
    tree_strings = f.readlines()
    index = 0
    tree_dict = dict()
    tree_index_dict = dict()
    for tree_str in tree_strings:
        tree_str = tree_str.split("'")[1]
        tree_dict[tree_str]=index
        tree_index_dict[index]=tree_str
        index += 1

    max_rank_moves = 0 # maximum number of rank moves among all shortest paths between all pairs of trees
    for i in range(0,len(tree_index_dict)):
        tree1 = read_from_cluster(tree_index_dict[i])
        for j in range(i+1,len(tree_index_dict)):
            tree2 = read_from_cluster(tree_index_dict[j])
            current_rank_moves = max_rank_move_shortest_path(tree1,tree2)
            if max_rank_moves < current_rank_moves:
                max_rank_moves = current_rank_moves
                print(tree_to_cluster_string(tree1), tree_to_cluster_string(tree2))
    return(max_rank_moves)


def test_bottom_up_hspr_approximation(num_leaves, hspr=1):
    if hspr == 1:
        d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves.npy')
        f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves.txt', 'r')
    else:
        d = np.load('SPR/distance_matrix_' + str(num_leaves) + '_leaves_hspr.npy')
        f = open('SPR/tree_dict_' + str(num_leaves) + '_leaves_hspr.txt', 'r')

    # Put all trees into a dict (note that indices are sorted increasingly in file)
    tree_strings = f.readlines()
    index = 0
    tree_dict = dict()
    tree_index_dict = dict()
    for tree_str in tree_strings:
        tree_str = tree_str.split("'")[1]
        tree_dict[tree_str]=index
        tree_index_dict[index]=tree_str
        index += 1
    
    # Most of the above is not necessary, as there seems to be a problem with dcd tr    

    differences = [] # array of differences between approximated and actual distances
    for i in range(0,len(d)):
        tree1 = read_from_cluster(tree_index_dict[i])
        for j in range(0,len(d)):
            tree2 = read_from_cluster(tree_index_dict[j])
            # print('start tree:', tree_to_cluster_string(tree1))
            # print('destination tree:', tree_to_cluster_string(tree2))
            # print(d[i][j], len(rankedspr_bfs(tree1, tree2, hspr=0))-1, rankedspr_path_bottom_up_hspr_dist(tree1,tree2))
            # print("computed path:")
            # path = rankedspr_path_bottom_up_hspr(tree1,tree2)
            # for i in range(0,path.num_trees):
            #     print(tree_to_cluster_string(path.trees[i]))
            # for some reason using the matrix d gives wrong results.
            # there might be something wrong with the computation of d??
            differences.append(rankedspr_path_bottom_up_hspr_dist(tree1,tree2)-(len(rankedspr_bfs(tree1, tree2, hspr=0))-1))
    print(differences)
