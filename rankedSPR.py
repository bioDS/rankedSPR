__author__ = 'Lena Collienne'
# Computing the rankedSPR graph to test algorithms for computing distances for trees on a small number of leaves
import sys
sys.path.append('treeOclock/')
sys.path.append('treeOclock/dct_parser/')

import ctypes
import math
import numpy as np
import random
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

    # Save adjacency matrix in file
    if hspr ==1:
        if not exists('SPR/adj_matrix_%s_leaves.npy' %num_leaves):
            np.save("SPR/adj_matrix_" + str(num_leaves) + "_leaves.npy", adj)
    else:
        if not exists('SPR/adj_matrix_%s_leaves_hspr.npy' %num_leaves):
            np.save("SPR/adj_matrix_" + str(num_leaves) + "_leaves_hspr.npy", adj)
    return(adj, tree_index)


def test_restricted_neighbourhood_search(num_leaves, num_tree_pairs):
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
        # print(tree_to_cluster_string(t_list.trees[i]))
        # print(tree_to_cluster_string(t_list.trees[i+1]))
        # print(rspr_distances[tree1_index][tree2_index], rankedspr_path_restricting_neighbourhood(t_list.trees[i],t_list.trees[i+1]))
        if (rspr_distances[tree1_index][tree2_index] == rankedspr_path_restricting_neighbourhood(t_list.trees[i],t_list.trees[i+1])):
            correct_distance += 1
    print('correct distance:', correct_distance, 'out of', num_tree_pairs)

test_restricted_neighbourhood_search(4,2)

# print(sum(rankedSPR_adjacency(4, hspr=0)[0]))
# print(sum(rankedSPR_adjacency(4)[0]))

# adj = np.load("SPR/adj_matrix_6_leaves.npy")
# print(adj.size)
# print(sum(sum(adj)))