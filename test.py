__author__ = 'Lena Collienne'
# Testing
import sys
sys.path.append('treeOclock/')
sys.path.append('treeOclock/dct_parser/')

import ctypes
import statistics as stat
import time
from treeOclock.dct_parser.tree_io import *
from treeOclock import *
from simulate_trees import *
from rankedSPR_seidel import *


t = "((1:1,2:1):2,(3:2,4:2):1);"
r = "(((1:1,2:1):1,3:2):1,4:3);"
s = "((3:1,4:1):2,(1:2,2:2):1);"

q1 = "((((1:1,2:1):1,4:2):1,5:3):1,3:4);"
q2 = "((((1:1,5:1):1,3:2):1,4:3):1,2:4);"

a = "(((1:1,2:1):1,3:2):3,((4:3,5:3):1,6:4):1);"
b = "(((4:1,5:1):1,6:2):3,((1:3,2:3):1,3:4):1);"

c = "(((1:1,2:1):2,3:3):2,((5:2,6:2):2,4:4):2);"
d = "(((5:1,6:1):2,4:3):2,((1:2,2:2):2,3:4):2);"

e = "((((1:1,2:1):1,3:2):1,4:3):3,((5:4,6:4):1,7:5):1)"
f = "((((1:3,2:3):1,3:4):1,4:5):1,((5:1,6:1):1,7:2):4)"


t1 = "((((3:1,4:1):1,5:2):2,6:4):1,(1:3,2:3):2);"
t2 = "((1:1,2:1):4,(((4:2,5:2):1,3:3):1,6:4):1);"

test = "(((2:1,3:1):2,1:3):1,(4:2,5:2):2);"


# ct = read_newick(q2)
# cs = read_newick(r)

# print(tree_to_cluster_string(ct))
# print(tree_to_cluster_string(del_leaf(ct, 1)))
# del_leaf(cs, 2)

distance_del_leaf(7,1,100,hspr=0)

# orbit_sizes(6,hspr=0)

# path = rankedspr_bfs(ct, cs, hspr=0)
# for i in range(0,len(path)):
#     print(path[i])

# ct = read_newick(t1, factor = 0)
# cr = read_newick(t2, factor = 0)
# print(ct, cs)
# print('distance(tree1, tree2):', len(rankedspr_bfs(ct,cs,hspr=0)))
# ct_neighbourhood = hspr_neighbourhood(ct)
# print('distance(neighbours of tree1, tree2):')
# for i in range(0, ct_neighbourhood.num_trees):
#     n_path = rankedspr_bfs(ct_neighbourhood.trees[i], cs, hspr=0)
#     if len(n_path)-1 == 2:
#         print("distance one neighbours:")
#         print(tree_to_cluster_string(ct_neighbourhood.trees[i]))
#         nested_neighbourhood = hspr_neighbourhood(ct_neighbourhood.trees[i])
#         print("distance two neighbours:")
#         for j in range(0, nested_neighbourhood.num_trees):
#             nn_path = rankedspr_bfs(nested_neighbourhood.trees[j], cs, hspr = 0)
#             if len(nn_path)-1 == 1:
#                 print(tree_to_cluster_string(nested_neighbourhood.trees[j]))


# compare_hspr_rspr(5,50)

# find_rank_moves(6,10)

# compare_hspr_rspr_uniform(6,100)

# print('symmetric cluster difference:', symmetric_cluster_diff(ct,cr))

# tree_list = sim_coal(5,100)
# for i in range(0, tree_list.num_trees - 2, 2):
#     print(i)
#     if len(rankedspr_bfs(tree_list.trees[i], tree_list.trees[i+1], 0))-1 != rankedspr_path_mrca_cluster_diff(tree_list.trees[i], tree_list.trees[i+1], 0).num_trees-1:
#         print("Not equal distance:")
#         print(tree_to_cluster_string(tree_list.trees[i]))
#         print(tree_to_cluster_string(tree_list.trees[i+1]))

# for tree in rankedspr_bfs(ct, cr, hspr=0):
#     print(tree)
# print('distance: ', len(rankedspr_bfs(ct, cr, 0))-1)

# path = rankedspr_path_mrca_cluster_diff(ct, cr, 0)
# print('approximated distance: ', path.num_trees-1)
# for i in range(0,path.num_trees):
#     print(tree_to_cluster_string(path.trees[i]))
#     print(symmetric_cluster_diff(path.trees[i],cr) + mrca_differences(path.trees[i],cr))

# rankedspr_path_restricting_neighbourhood(ct,cr)

# rankedspr_seidel(4,0)
# orbit_size_identity_caterpillar(4)

# a = "((((3:1,4:1):1,5:2):2,6:4):1,(1:3,2:3):2);"
# b = "((1:1,2:1):4,(((4:2,5:2):1,3:3):1,6:4):1);"

# c = "((1:1,2:1):5,((((3:2,4:2):1,5:3):1,6:4):1,7:5):1)"
# d = "((1:1,2:1):5,((3:2,7:2):3,((4:3,5:3):1,6:4):1):1);"

# t1 = "(((2:1,3:1):2,(1:2,4:2):1):1,5:4);"
# t2 = "((((1:1,5:1):1,4:2):1,3:3):1,2:4);"

# t3 = "(((1:1,5:1):2,(3:2,4:2):1):1,2:4);"
# t4 = "((((1:1,5:1):1,4:2):1,2:3):1,3:4);"

# t5 = "(((2:1,5:1):1,3:2):2,(1:3,4:3):1);"
# t6 = "(((1:1,5:1):2,3:3):1,(2:2,4:2):2);"


# tree1 = read_newick(t5, factor=0)
# tree2 = read_newick(t6, factor=0)

# print("Top-down approach:")
# print("approx distance:", rankedspr_path_top_down_symm_diff(tree1, tree2))

# print("Actual shortest path:")
# path = rankedspr_bfs(tree1, tree2)
# print('distance:', len(path)-1)
# print('path:')
# for tree in path:
#     print(tree)

# print(symm_cluster_diff(tree1, tree2, 2*tree1.num_leaves-3))

# test_top_down_neighbourhood_search(5,100)

# n=7
# m=1000

#caterpillar_trees = sim_cat(n,m)
#identity_caterpillar = identity_caterpillar(n)

# for i in range(0,m):
#     print(i)
#     path = rankedspr_bfs(identity_caterpillar, caterpillar_trees.trees[i])
#     if len(path) >7:
#         print('path:')
#         for tree in path:
#             print(tree)

# for n in range(4,7):
#     print('n =',n)
#     # rankedspr_seidel(n, hspr = 1)
#     caterpillar_diameter_trees(n)

# n=7
# m =100
#
# # PROBLEM: only a very small number of trees actually has diameter distance from each other.
# # For 7 leaves we are very unlikely to get those by chance.
#for i in range(0,m):
#    sim_trees = sim_coal(n, 2) #simulate 2 ranked trees with n leaves, repeat this m times
#    # compare pairwise distances between trees the two tree
#    path = rankedspr_bfs(sim_trees.trees[0], sim_trees.trees[1])
#    print(i, 'distance:', len(path)-1)
#    if len(path) - 1 == 7:
#        print("path:")
#        for tree in path:
#            print(tree)


######## TESTING ######## (moved here from rankedSPR.py)

# coal_pw_spr_dist(5, 10000, hspr = 0)
# coal_pw_spr_dist(5, 10000, hspr = 1)


# test_restricted_neighbourhood_search(5,1000, hspr = 1)
# test_top_down_neighbourhood_search(5,10)

# print(sum(rankedSPR_adjacency(4, hspr=0)[0]))
# print(sum(rankedSPR_adjacency(4)[0]))


# t1 = "(((1:1,2:1):2,(3:2,4:2):1):1,5:4);"
# t2 = "((((3:1,5:1):1,2:2):1,1:3):1,4:4);"


# identity_caterpillar_7 = "((((((1:1,2:1):1,3:2):1,4:3):1,5:4):1,6:5):1,7:6);"
# reverse_identity_caterpillar_7 = "((((((6:1,7:1):1,5:2):1,4:3):1,3:4):1,2:5):1,1:6);"

# identity_caterpillar_6 = "(((((1:1,2:1):1,3:2):1,4:3):1,5:4):1,6:5);"
# reverse_identity_caterpillar_6 = "(((((5:1,6:1):1,4:2):1,3:3):1,2:4):1,1:5);"

# identity_caterpillar_5 = "((((1:1,2:1):1,3:2):1,4:3):1,5:4);"
# reverse_identity_caterpillar_5 = "((((4:1,5:1):1,3:2):1,2:3):1,1:4);"

# identity_caterpillar_4 = "(((1:1,2:1):1,3:2):1,4:3);"
# reverse_identity_caterpillar_4 = "(((3:1,4:1):1,2:2):1,1:3);"

# tree1 = read_newick(identity_caterpillar_4, factor = 0)
# tree2 = read_newick(reverse_identity_caterpillar_4, factor = 0)

# path = rankedspr_bfs(tree1, tree2)
# print("distance:", len(path)-1)

# for tree in path:
#     print(tree)

# print("rankedSPR neighbours:")
# spr_neighbours = spr_neighbourhood(tree1)
# for i in range(0,spr_neighbours.num_trees):
#     print(tree_to_cluster_string(spr_neighbours.trees[i]))


# Compute distance matrices in HSPR and RSPR for small number of leaves and print number of trees with diameter distance
# for n in range(4,7):
#     rankedspr_seidel(n,0)
#     rankedspr_seidel(n,1)

#     print("RSPR for " + str(n) + " leaves:")
#     dist = np.load("SPR/distance_matrix_" + str(n) + "_leaves.npy")
#     print(len(np.where(dist == np.amax(dist))[0])/2)

#     print("HSPR for " + str(n) + " leaves:")
#     dist = np.load("SPR/distance_matrix_" + str(n) + "_leaves_hspr.npy")
#     print(len(np.where(dist == np.amax(dist))[0])/2)

# adj = np.load("SPR/adj_matrix_5_leaves.npy")
# tree1_1nh = set()
# tree2_2nh = set()
# print("adjacency tree1:")
# for i in range(0,180):
#     if(adj[7][i] == 1):
#         print(i)
#         tree1_1nh.add(i)
# print("adjacency tree2:")
# for i in range(0,180):
#     if(adj[154][i] == 1):
#         for j in range(0,180):
#             if(adj[i][j]==1):
#                 tree2_2nh.add(j)
# print(tree1_1nh)
# print(tree2_2nh)
# print(tree1_1nh.intersection(tree2_2nh))
# print(adj[7])
# print(adj[154])
# print(adj.size)
# print(sum(sum(adj)))
