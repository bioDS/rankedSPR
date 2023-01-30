// Implementation of algorithm computing HSPR path as in Theorem 1 of ranked SPR paper

#include "spr_path.h"

// Create a path by using a bottom-up approach in RSPR/HSPR, only using HSPR
// moves -- greedily create dest_tree from bottom to top
Tree_Array rankedspr_path_bottom_up_hspr(Tree* start_tree, Tree* dest_tree) {
    long num_leaves = start_tree->num_leaves;

    // Initialise output path
    Tree_Array path = get_empty_tree_array((num_leaves - 1) * (num_leaves - 2) +
                     1, num_leaves);

    // Deep copy start tree to get new tree to be added to path iteratively
    Tree* current_tree = new_tree_copy(start_tree);

    // Add the first tree to output path
    for (long i = 0; i < 2 * num_leaves - 1; i++) {
        path.trees[0].node_array[i] = current_tree->node_array[i];
    }

    long index = 1;  // index of the currently last tree on output path

    for (long i = num_leaves; i < 2 * num_leaves - 1; i++) {
        long current_child1 = current_tree->node_array[i].children[0];
        long current_child2 = current_tree->node_array[i].children[1];
        if ((current_child1 == dest_tree->node_array[i].children[0] &&
             current_child2 == dest_tree->node_array[i].children[1]) ||
            (current_child2 == dest_tree->node_array[i].children[0] &&
             current_child1 == dest_tree->node_array[i].children[1])) {
            // don't do anything, proceed to next iteration
        } else if (current_child1 == dest_tree->node_array[i].children[0]) {
            spr_move(current_tree, i, dest_tree->node_array[i].children[1], 0);
            // add current_tree to path
            for (long j = 0; j < 2 * num_leaves - 1; j++) {
                path.trees[index].node_array[j] = current_tree->node_array[j];
            }
            index++;
        } else if (current_child1 == dest_tree->node_array[i].children[1]) {
            spr_move(current_tree, i, dest_tree->node_array[i].children[0], 0);
            // add current_tree to path
            for (long j = 0; j < 2 * num_leaves - 1; j++) {
                path.trees[index].node_array[j] = current_tree->node_array[j];
            }
            index++;
        } else if (current_child2 == dest_tree->node_array[i].children[0]) {
            spr_move(current_tree, i, dest_tree->node_array[i].children[1], 1);
            // add current_tree to path
            for (long j = 0; j < 2 * num_leaves - 1; j++) {
                path.trees[index].node_array[j] = current_tree->node_array[j];
            }
            index++;
        } else if (current_child2 == dest_tree->node_array[i].children[1]) {
            spr_move(current_tree, i, dest_tree->node_array[i].children[0], 1);
            // add current_tree to path
            for (long j = 0; j < 2 * num_leaves - 1; j++) {
                path.trees[index].node_array[j] = current_tree->node_array[j];
            }
            index++;
        } else {  // choose a smaller child (lower index) of current node i to
                  // move << This choice of child influences the length of the
                  // output path!
            int child_index_dest_tree = 0;
            if (dest_tree->node_array[i].children[0] >
                dest_tree->node_array[i].children[1]) {
                child_index_dest_tree = 1;
            }
            int child_index_current_tree = 0;
            if (current_child1 > current_child2) {
                child_index_current_tree = 1;
            }
            spr_move(current_tree, i,
                     dest_tree->node_array[i].children[child_index_dest_tree],
                     child_index_current_tree);
            // add current_tree to path
            for (long j = 0; j < 2 * num_leaves - 1; j++) {
                path.trees[index].node_array[j] = current_tree->node_array[j];
            }
            index++;
            // find the index of the child that we want to move now (the one
            // that coincides with the one in dest_tree)
            int child_index = 0;
            if (current_tree->node_array[i].children[1] ==
                dest_tree->node_array[i].children[0]) {
                child_index = 1;
            }
            // now move the child that has correct parent to the other one of
            // dest_tree
            spr_move(
                current_tree, i,
                dest_tree->node_array[i].children[1 - child_index_dest_tree],
                child_index);
            // add current_tree to path
            for (long j = 0; j < 2 * num_leaves - 1; j++) {
                path.trees[index].node_array[j] = current_tree->node_array[j];
            }
            index++;
        }
    }
    path.num_trees = index;
    free(current_tree);
    return (path);
}