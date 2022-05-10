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


n=7
m=100

caterpillar_trees = sim_cat(n,m)

for i in range(0,m,2):
    print(i)
    path = rankedspr_bfs(caterpillar_trees.trees[i], caterpillar_trees.trees[i+1])
    if len(path) >7:
        print('path:')
        for tree in path:
            print(tree)


# n=5
# m=100
# sim_trees = sim_coal(n, 100, rf = False) #simulate m ranked trees with n leaves

# for i in range(0,99,2):
#     # compare pairwise distances between trees at positions i and i+1 in sim_trees
#     path = rankedspr_path_mrca_cluster_diff(sim_trees.trees[i], sim_trees.trees[i+1])
#     if path.num_trees - 1 >4:
#         for i in range(0,path.num_trees):
#             print(tree_to_cluster_string(path.trees[i]))