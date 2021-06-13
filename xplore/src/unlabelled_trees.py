__author__ = 'Lena Collienne'
# Unlabelled RNNI

from os import unlink
import copy
import os.path
import sys
from ete3.treeview.main import NODE_STYLE_DEFAULT
# sys.path.append('../..')

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
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
    u_tree = copy.deepcopy(unlabelled_tree) # Copy unlabelled tree, as we will pop elements of the sets in the list
    num_leaves = len(u_tree) + 2
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
        child_1 = u_tree[i].pop()
        if (len(u_tree[i]) > 0): #the set might be empty if there was only one zero in there
            child_2 = u_tree[i].pop()
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


# Return a labelled tree such that leaf labels increase with increasing rank
def label_tree_increasingly(unlabelled_tree):
    u_tree = copy.deepcopy(unlabelled_tree) # Copy unlabelled tree, as we will pop elements of the sets in the list
    num_leaves = len(u_tree) + 2
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
        child_1 = u_tree[i].pop()
        if (len(u_tree[i]) > 0): #the set might be empty if there was only one zero in there
            child_2 = u_tree[i].pop()
        else:
            child_2 = 0
        # SOMETHING HERE IS QUITE BROKEN!!!!
        if child_1 != 0:
            node_list[num_leaves+i+1].children[0] = num_leaves-1+child_1
            node_list[num_leaves+child_1-1].parent = num_leaves+i+1
        else:
            leaf = leaves.pop() # choose last leaf label as child (highest possible leaf label)
            node_list[num_leaves+i+1].children[0] = leaf
            node_list[leaf].parent = num_leaves+i+1
        if child_2 != 0:
            node_list[num_leaves+i+1].children[1] = num_leaves-1+child_2
            node_list[num_leaves+child_2-1].parent = num_leaves+i+1
        else:
            leaf = leaves.pop() # choose last leaf label as child (highest possible leaf label)
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


def unlabelled_dist(Tlist, Rlist):
    dist = 0
    for i in range(0,len(Tlist)):
        if Tlist[i] != {0} or Rlist[i] != {0}:
            dist += 2-len(Tlist[i].intersection(Rlist[i]))
    return(dist)


def approx_unlabelled_RNNI_dist(t1, t2, N):
    # Approximate the unlabelled distance between unlabelled trees t1 and t2 (given in unlabelled tree list representation) by randomly labelling t1 (fixed) and then trying N different labellings for t2
    tree1 = unlabelled_to_labelled_tree(t1)
    dist_list = []
    for i in range(0, N):
        tree2 = unlabelled_to_labelled_tree(t2)
        dist_list.append(findpath_distance(tree1, tree2))
    return(min(dist_list))


def compare_arbitrary_to_increasing_labelling_URNNI(n,m,l,ldist = False):
    # Simulate m trees on n leaves (coalescent, then delete labels) and compare the URNNI dist proxies resulting from labelling increasingly with rank and approx_labelling from l arbitrary labellings
    incr_label_dist = []
    arb_label_dist = []
    if ldist == True:
        list_dist = []
    diff = []
    for i in range(0,m):
        tree_list = sim_coal(n,2)
        t1 = labelled_to_unlabelled_tree(tree_list.trees[0])
        t2= labelled_to_unlabelled_tree(tree_list.trees[1])
        incr_label_dist.append(findpath_distance(label_tree_increasingly(t1), label_tree_increasingly(t2)))
        arb_label_dist.append(approx_unlabelled_RNNI_dist(t1,t2,l))
        if ldist == True:
            list_dist.append(unlabelled_dist(t1,t2))
        diff.append(incr_label_dist[i] - arb_label_dist[i])
    # Plot distances or differences btw distances
    if ldist == True:
        d = pd.DataFrame(data = list(zip(incr_label_dist, arb_label_dist, list_dist)), columns = ["increasing labelling", "arbitrary labelling", "list_dist"])
    else:
        d = pd.DataFrame(data = diff)
    sns.scatterplot(data=d, legend = True)
    # plt.tight_layout()
    plt.legend(bbox_to_anchor=(0.6, 0.3), loc='upper left', borderaxespad=0)
    if ldist == True:
        plt.savefig("unlabelled_RNNI_plots/compare_incr_arbr_labellings_list_dist/" + str(n) + "_leaves_" + str(m) + "_repeats_" + str(l) +  "_arbitrary_labellings.pdf")
    else:
        plt.savefig("unlabelled_RNNI_plots/compare_labellings/" + str(n) + "_leaves_" + str(m) + "_repeats_" + str(l) +  "_arbitrary_labellings.pdf")
    plt.show()


def compare_increasing_labellings_list_dist(n,m):
    # Simulate m trees on n leaves (coalescent, then delete labels) and compare the URNNI dist proxy resulting from labelling increasingly with rank to list dist for unlabelled trees
    incr_label_dist = []
    list_dist = []
    diff = []
    for i in range(0,m):
        tree_list = sim_coal(n,2)
        t1 = labelled_to_unlabelled_tree(tree_list.trees[0])
        t2= labelled_to_unlabelled_tree(tree_list.trees[1])
        incr_label_dist.append(findpath_distance(label_tree_increasingly(t1), label_tree_increasingly(t2)))
        list_dist.append(unlabelled_dist(t1,t2))
        diff.append(incr_label_dist[i] - list_dist[i])
    # Plot 
    d = pd.DataFrame(data = list(zip(incr_label_dist, list_dist)), columns = ["increasing labelling", "list labelling"])
    # d = pd.DataFrame(data = diff)
    sns.scatterplot(data=d, legend = True)
    plt.tight_layout()
    # plt.savefig("unlabelled_RNNI_plots/compare_incr_labelling_list_dist/plot_both" + str(n) + "_leaves_" + str(m) + "_repeats.pdf")
    # plt.savefig("unlabelled_RNNI_plots/compare_incr_labelling_list_dist/" + str(n) + "_leaves_" + str(m) + "_repeats.pdf")
    plt.show()



if __name__ == '__main__':

    compare_arbitrary_to_increasing_labelling_URNNI(100,100,1000, ldist = True)
    # compare_increasing_labellings_list_dist(100,100)

    # labelled_tree = sim_coal(20,1).trees[0]
    # t1 = [{0,1}, {0,0}, {2,3}]
    # t2 = [{0,0}, {1,0}, {2,3}]
    # print(findpath_distance(label_tree_increasingly(t1), label_tree_increasingly(t2)))
    # print(unlabelled_dist(t1, t2))
    # # labelled_to_unlabelled_tree(c_tree)
    # t = unlabelled_to_labelled_tree(t1)
    # utree = labelled_to_unlabelled_tree(t)
    # print(t1, t2)
    # print(unlabelled_dist(utree,utree))
    # for i in range(0,20):
    #     apprx_RNNI = approx_unlabelled_RNNI_dist(t1,t2,1000)
    #     u_dist = unlabelled_dist(t1,t2)
    #     print(u_dist/apprx_RNNI)
