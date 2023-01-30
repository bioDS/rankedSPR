__author__ = 'Lena Collienne'
# Simulating Trees under coalescent and birth-death(yule) model.
# Trees genereted in both DCT and ete3 format.

import sys


from treeOclock.tree_parser.tree_io import tree_to_cluster_string

import random
from itertools import permutations
import math

from tree_functions import *

def sim_coal(num_leaves, num_trees):
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
            if j >= num_leaves:
                node_list[j].time = j-num_leaves+1
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
    output_tree_list = TREE_ARRAY(trees, num_trees)
    return(output_tree_list)


def del_leaf(tree, leaf):
    # deletes leaf at position 'leaf' from tree

    # Create empty Node list for output tree on n-1 leaves
    num_leaves = tree.num_leaves - 1
    num_nodes = 2*num_leaves - 1 #new tree will have n-1 leaves and n-2 internal nodes
    node_list = (NODE * num_nodes)()
    empty_children = (c_long * 2)()
    empty_children[0] = -1
    empty_children[1] = -1
    # initialise node list:
    for i in range(0, num_nodes):
        node_list[i] = NODE(-1, empty_children, 0)
        if i >= num_leaves:
            node_list[i].time = i-num_leaves+1
    for i in range(0, num_nodes-1):
        # i is the index in node_list, k is index in tree
        if i < leaf:
            k = i
        elif i < tree.node_array[leaf].parent-1:
            k = i+1
        else:
            k = i+2
        # print("i:", i, "k:", k)
        if k == tree.node_array[leaf].parent: # update parent of the node that used to be parent of the deleted leaf
            node_list[i].parent = tree.node_array[tree.node_array[leaf].parent].parent - 2
        elif tree.node_array[k].parent == tree.node_array[leaf].parent: # update parent for sibling of deleted leaf
            node_list[i].parent = tree.node_array[tree.node_array[leaf].parent].parent - 2
        elif tree.node_array[k].parent > tree.node_array[leaf].parent:
            node_list[i].parent = tree.node_array[k].parent-2 # subtract two for one leaf and one internal node less
        elif tree.node_array[k].parent == -1:
            node_list[i].parent = -1
        else:
            node_list[i].parent = tree.node_array[k].parent-1 # subtract one for one leaf less

    # now update children
    for i in range(num_leaves, num_nodes):
        if i < tree.node_array[leaf].parent-1:
            k = i+1
        else:
            k = i+2
        for j in [0,1]:
            p = tree.node_array[leaf].parent
            if tree.node_array[k].children[j] == p:
                # replace the child that is 'leaf' by the sibling of 'leaf' in the original tree
                if tree.node_array[p].children[0] == leaf:
                    if tree.node_array[p].children[1] < leaf:
                        node_list[i].children[j] = tree.node_array[p].children[1]
                    else:
                        node_list[i].children[j] = tree.node_array[p].children[1]-1
                else:
                    if tree.node_array[p].children[0] < leaf:
                        node_list[i].children[j] = tree.node_array[p].children[0]
                    else:
                        node_list[i].children[j] = tree.node_array[p].children[0]-1
            elif tree.node_array[k].children[j] < leaf:
                node_list[i].children[j] = tree.node_array[k].children[j]
            elif tree.node_array[k].children[j] < p:
                node_list[i].children[j] = tree.node_array[k].children[j]-1
            else:
                node_list[i].children[j] = tree.node_array[k].children[j]-2

    output = TREE(node_list, num_leaves)
    return(output)


def identity_caterpillar(num_leaves):
    # Create empty node list:
    num_nodes = 2*num_leaves-1
    node_list = (NODE * num_nodes)()
    empty_children = (c_long * 2)()
    empty_children[0] = -1
    empty_children[1] = -1
    for j in range(0, num_nodes):
        node_list[j] = NODE(-1, empty_children, 0)
        if j >= num_leaves:
            node_list[j].time = j-num_leaves+1
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