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
from rankedSPR import *


t = "((1:1,2:1):2,(3:2,4:2):1);"
r = "(((1:1,2:1):1,3:2):1,4:3);"
s = "(((3:1,4:1):1,2:2):1,1:3);"

q1 = "((((1:1,2:1):1,4:2):1,5:3):1,3:4);"
q2 = "((((4:1,5:1):1,3:2):1,2:3):1,1:4);"

# ct = read_newick(q1, factor = 0)
# cr = read_newick(q2, factor = 0)

# print('symmetric cluster difference:', symmetric_cluster_diff(ct,cr))

# path = rankedspr_path_mrca_cluster_diff(ct, cr)
# print('distance: ', path.num_trees-1)
# for i in range(0,path.num_trees):
#     print(tree_to_cluster_string(path.trees[i]))
#     print(symmetric_cluster_diff(path.trees[i],cr) + mrca_differences(path.trees[i],cr))

# rankedspr_path_restricting_neighbourhood(ct,cr)


# n=7
# m=1000

# caterpillar_trees = sim_cat(n,m)
# identity_caterpillar = identity_caterpillar(n)

# for i in range(0,m):
#     print(i)
#     path = rankedspr_bfs(identity_caterpillar, caterpillar_trees.trees[i])
#     if len(path) >7:
#         print('path:')
#         for tree in path:
#             print(tree)

# d = np.load('SPR/distance_matrix_7_leaves.npy')
# print(np.count_nonzero(d==7))

n=7
m=100

PROBLEM: only a very small number of trees actually has diameter distance from each other.
For 7 leaves we are very unlikely to get those by chance.
for i in range(0,m):
    sim_trees = sim_coal(n, 2) #simulate 2 ranked trees with n leaves, repeat this m times
    # compare pairwise distances between trees the two tree
    path = rankedspr_bfs(sim_trees.trees[0], sim_trees.trees[1])
    print(i, 'distance:', len(path)-1)
    if len(path) - 1 == 7:
        print("path:")
        for tree in path:
            print(tree)
