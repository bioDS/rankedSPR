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

def sim_cat(num_leaves, num_trees):
    # simulate num_trees caterpillar trees (random at uniform)
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
        # First simulate cherry
        [c1,c2] = random.sample(current_leaves, k=2)
        current_leaves.remove(c1)
        current_leaves.remove(c2)
        node_list[c1-1].parent = num_leaves
        node_list[c2-1].parent = num_leaves
        node_list[num_leaves].children[0] = c1-1
        node_list[num_leaves].children[1] = c2-1
        # Add leaves above cherry to caterpillar tree
        for node in range(2,num_leaves):
            leaf = random.choice(current_leaves) # Choose next leaf to add random at uniform
            current_leaves.remove(leaf)
            # Add relations in node_list
            node_list[leaf-1].parent = num_leaves + node - 1
            node_list[num_leaves+node-2].parent = num_leaves+node-1
            node_list[num_leaves+node-1].children[0] = leaf-1
            node_list[num_leaves+node-1].children[1] = num_leaves+node-2
        # for l in range(0, num_nodes):
        #     print(l, node_list[l].children[0], node_list[l].children[1], node_list[l].parent)
        current_tree = TREE(node_list, num_leaves)
        trees[i] = current_tree
    output_tree_list = TREE_LIST(trees, num_trees)
    return(output_tree_list)

def balanced_tree_16_leaves():
    num_leaves = 16
    num_nodes = 2 * num_leaves - 1
    # Create empty Node list
    node_list = (NODE * num_nodes)()
    empty_children = (c_long * 2)()
    empty_children[0] = -1
    empty_children[1] = -1
    for j in range(num_leaves, num_nodes):
        node_list[j] = NODE(-1, empty_children, 0)
    # set all parents of leaves
    for i in range(num_leaves, int(num_leaves + num_leaves/2), 2):
        node_list[i].children[0] = i-num_leaves
        node_list[i].children[1] = i-num_leaves+1
        node_list[i-num_leaves].parent = i
        node_list[i-num_leaves+1].parent = i
    for i in range(num_leaves + 1, int(num_leaves + num_leaves/2), 2):
        node_list[i].children[0] = int(i-num_leaves/2)-1
        node_list[i].children[1] = int(i-num_leaves/2)
        node_list[int(i-num_leaves/2)-1].parent = i
        node_list[int(i-num_leaves/2)].parent = i
    # Set parents of internal nodes
    # Everything involving root
    node_list[num_nodes-1].children[0] = num_nodes-2
    node_list[num_nodes-1].children[1] = num_nodes-3
    node_list[num_nodes-2].parent = num_nodes-1
    node_list[num_nodes-3].parent = num_nodes-1

    # Children of root
    node_list[num_nodes-2].children[0] = num_nodes-6
    node_list[num_nodes-2].children[1] = num_nodes-4
    node_list[num_nodes-3].children[0] = num_nodes-7
    node_list[num_nodes-3].children[1] = num_nodes-5
    node_list[num_nodes-6].parent = num_nodes-2
    node_list[num_nodes-4].parent = num_nodes-2
    node_list[num_nodes-7].parent = num_nodes-3
    node_list[num_nodes-5].parent = num_nodes-3

    # Children of children of root
    node_list[num_nodes-6].children[0] = num_leaves + 1
    node_list[num_nodes-6].children[1] = num_leaves + 3
    node_list[num_nodes-4].children[0] = num_leaves + 5
    node_list[num_nodes-4].children[1] = num_leaves + 7
    node_list[num_nodes-7].children[0] = num_leaves
    node_list[num_nodes-7].children[1] = num_leaves + 2
    node_list[num_nodes-5].children[0] = num_leaves + 4
    node_list[num_nodes-5].children[1] = num_leaves + 6
    node_list[num_leaves+1].parent = num_nodes-6
    node_list[num_leaves+3].parent = num_nodes-6
    node_list[num_leaves+5].parent = num_nodes-4
    node_list[num_leaves+7].parent = num_nodes-4
    node_list[num_leaves].parent = num_nodes-7
    node_list[num_leaves+2].parent = num_nodes-7
    node_list[num_leaves+4].parent = num_nodes-5
    node_list[num_leaves+6].parent = num_nodes-5
    output_tree = TREE(node_list, num_leaves)
    return(output_tree)