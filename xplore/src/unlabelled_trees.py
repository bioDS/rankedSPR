__author__ = 'Lena Collienne'
# Unlabelled RNNI

from os import unlink
import os.path
import sys
from ete3.treeview.main import NODE_STYLE_DEFAULT
# sys.path.append('../..')

import numpy as np
import random

from simulate_trees import *

# Return the underlying unlabelled tree for a given labelled tree
# Note: If a set only contains one 0, the corresponding internal node has two leaves as children (because we use sets)
def labelled_to_unlabelled_tree(tree):
    num_leaves = tree.num_leaves
    # Initialise unlabelled_tree (output) as list of empty sets (these will contain the two children for each internal node in the end)
    unlabelled_tree = []
    for i in range(1,num_leaves-1):
        unlabelled_tree.append(set())
    for i in range(2*num_leaves-2, num_leaves, -1):
        # top-down: Loop from root to leaves
        child_1 = tree.tree[i].children[0]
        child_2 = tree.tree[i].children[1]
        # print(i, child_1, child_2)
        if child_1 > num_leaves-1: # Check if child_1 of current node is an internal node
            unlabelled_tree[i-num_leaves-1].add(child_1-num_leaves+1)
        else:
            unlabelled_tree[i-num_leaves-1].add(0)
        if child_2 > num_leaves-1: # Check if child_2 of current node is an internal node
            # print(i, child_2)
            unlabelled_tree[i-num_leaves-1].add(child_2-num_leaves+1)
            # print(i, child_2, i-num_leaves-1)
            # print(unlabelled_tree[i-num_leaves-1])
        else:
            unlabelled_tree[i-num_leaves-1].add(0)
    return(unlabelled_tree)

# Return a labelled tree (randomly assigned leaf labels)
def unlabelled_to_labelled_tree(unlabelled_tree):
    u_tree = unlabelled_tree # Copy unlabelled tree, as we will pop elements of the sets in the list
    num_leaves = len(unlabelled_tree) + 2
    # Initialise output tree:
    num_nodes = 2 * num_leaves - 1
    # Create empty Node list
    node_list = (NODE * num_nodes)()
    empty_children = (c_long * 2)()
    empty_children[0] = -1
    empty_children[1] = -1
    for j in range(0, num_nodes):
        node_list[j] = NODE(-1, empty_children, 0)
    # Set of children, which will be assigned randomly at places where unlabelled tree representation has zeros:
    leaves = list()
    for i in range(0, num_leaves):
        leaves.append(i)
    for i in range(num_leaves-3, -1, -1):
        # Top-down approach: fill tree from root to leaves
        child_1 = unlabelled_tree[i].pop()
        if (len(unlabelled_tree[i]) > 0): #the set might be empty if there was only one zero in there
            child_2 = unlabelled_tree[i].pop()
        else:
            child_2 = 0
        # SOMETHING HERE IS QUITE BROKEN!!!!
        if child_1 != 0:
            node_list[num_leaves+i+1].children[0] = num_leaves-1+child_1
            node_list[num_leaves+child_1-1].parent = num_leaves+i+1
        else:
            r = random.randint(0,len(leaves)-1)
            leaf = leaves[r]
            leaves.pop(r) # choose an arbitrary leaf label as child
            node_list[num_leaves+i+1].children[0] = leaf
            node_list[leaf].parent = num_leaves+i+1
        if child_2 != 0:
            node_list[num_leaves+i+1].children[1] = num_leaves-1+child_2
            node_list[num_leaves+child_2-1].parent = num_leaves+i+1
        else:
            r = random.randint(0,len(leaves)-1)
            leaf = leaves[r]
            leaves.pop(r) # choose an arbitrary leaf label as child
            node_list[num_leaves+i+1].children[1] = leaf
            node_list[leaf].parent = num_leaves+i+1
    # Fill cherry of rank 1:
    child_1 = leaves.pop()
    child_2 = leaves.pop()
    node_list[num_leaves].children[0]=child_1
    node_list[num_leaves].children[1]=child_2
    node_list[child_1].parent = num_leaves
    node_list[child_2].parent = num_leaves

    # Check if the tree is read in correctly:
    # for i in range(len(node_list)-1,-1,-1):
    #     print(node_list[i].children[0], node_list[i].children[1], node_list[i].parent)
    
    output = TREE(node_list, num_leaves)
    return(output)

labelled_tree = sim_coal(20,1).trees[0]
unlabelled_tree = [{0,1}, {0,0}, {2,3}]
# labelled_to_unlabelled_tree(c_tree)
t = unlabelled_to_labelled_tree(unlabelled_tree)
print(labelled_to_unlabelled_tree(t))