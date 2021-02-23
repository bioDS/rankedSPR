__author__ = 'Lena Collienne'
# Simulating Trees under coalescent and birth-death(yule) model.
# Trees genereted in both DCT and ete3 format.

import os.path
import sys
sys.path.append('../..')

from ete3 import Tree
import numpy as np
import random

from tree_structs import *
from dct_parser.tree_io import *

def sim_coal(num_leaves, num_trees, rf = False):
    trees = (TREE * num_trees)()
    num_nodes = 2 * num_leaves - 1
    for i in range(0, num_trees):
        # Create empty Node list
        node_list = (NODE * num_nodes)()
        empty_children = (c_long * 2)()
        empty_children[0] = -1
        empty_children[1] = -1
        for j in range(0, num_nodes):
            node_list[j] = NODE(-1, empty_children, 0)
        # Simulate coalescence events
        n = num_leaves # n decreases in every step of simulation until all lineages coalesce
        current_leaves = [] # list including all current leaves. Starts with list of leaves, in each iteration two elements are replaced by one (which is attached as last element of the list)
        for l in range(1,num_leaves+1):
            current_leaves.append(l)
        for j in range(0,num_leaves-1):
            [n1,n2] = random.sample(current_leaves, k=2)
            # leaves coalescing in internal node of rank j+1 (index n+j+1 in node_list)
            current_leaves.remove(n1)
            current_leaves.remove(n2)
            current_leaves.append(int(num_leaves+j+1)) #add new internal node to leaf set
            # Add new relations to node_list:
            node_list[n1-1].parent = num_leaves + j
            node_list[n2-1].parent = num_leaves + j
            node_list[num_leaves+j].children[0] = n1-1
            node_list[num_leaves+j].children[1] = n2-1
        # for l in range(0, num_nodes):
        #     print(l, node_list[l].children[0], node_list[l].children[1], node_list[l].parent)
        current_tree = TREE(node_list, num_leaves)
        trees[i] = current_tree
    output_tree_list = TREE_LIST(trees, num_trees)
    return(output_tree_list)

def identity_caterpillar(num_leaves):
    # Create empty node list:
    num_nodes = 2*num_leaves-1
    node_list = (NODE * num_nodes)()
    empty_children = (c_long * 2)()
    empty_children[0] = -1
    empty_children[1] = -1
    for j in range(0, num_nodes):
        node_list[j] = NODE(-1, empty_children, 0)
    # First set cherry, remaining leaves get parents in loop
    node_list[0].parent = num_leaves
    node_list[1].parent = num_leaves
    node_list[num_leaves].children[0] = 0
    node_list[num_leaves].children[1] = 1
    for leaf in range(2,num_leaves):
        node_list[leaf].parent = num_leaves + leaf - 1
        node_list[num_leaves + leaf - 1].children[0] = leaf
        node_list[num_leaves + leaf - 1].children[1] = num_leaves + leaf - 2
        node_list[num_leaves + leaf - 2].parent = num_leaves + leaf - 1
    # for i in range(0, num_nodes):
    #     print(node_list[i].children[0], node_list[i].children[1], node_list[i].parent)
    output_tree = TREE(node_list, num_leaves)
    return(output_tree)
