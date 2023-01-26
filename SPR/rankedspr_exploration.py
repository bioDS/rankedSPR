__author__ = 'Lena Collienne'
# Exploring some properties of the RSPR and HSPR graph

from rankedSPR_seidel import *
from SPR.rankedspr_adjacency import *
from plots import *
from spr_playground_functions import *
import gc


def add_unique_topology_pair(tree_pairs, tree1, tree2):
    # only add tree1 and tree2 to tree_pairs if permuting leaves of tree1 and tree2 in the same way does not give a tree pair in tree_pairs.
    tp_exists = False
    for tp in tree_pairs:
        tree3_str = min(tp)
        tree4_str = max(tp)
        tree3 = read_from_cluster(tree3_str)
        tree4 = read_from_cluster(tree4_str)
        if tree_pairs_permutations(tree1, tree2, tree3, tree4):
            tp_exists = True
            break
    if tp_exists == False:
        tree1_str = tree_to_cluster_string(tree1)
        tree2_str = tree_to_cluster_string(tree2)
        tree_pairs.append(set([tree1_str, tree2_str]))
    return tree_pairs


# Find the ranks on which moves are performed on shortest path resulting from BFS
def bfs_path_rank_sequence(tree1, tree2):
    path = rankedspr_bfs(tree1, tree2, hspr=True)
    rank_list = []
    rank_count = []
    for i in range(0, len(path) - 1):
        # for each move, find the lowest rank for which induced cluster changes -- this is the rank on which the HSPR move happened
        path[i] = str(path[i])
        path[i + 1] = str(path[i + 1])
        rank = 0
        for j in range(0, len(path[i]) - 1):
            if path[i][j] == '{':
                rank += 1
            if path[i][j] != path[i + 1][j]:
                rank_list.append(rank)
                break
    for i in range(1, tree1.num_leaves - 1):
        rank_count.append(rank_list.count(i))
    if max(rank_count) > 3:
        print("There is a rank for which more than one move is needed. The corresponding path is:")
        for tree in path:
            print(tree)
    return (rank_count)


def check_HSPR_moves_per_rank(num_leaves, num_tree_pairs):
    # simulate num_tree_pairs trees and check how moves are distributed across ranks in the trees on shortest path computed by bfs
    for i in range(0, num_tree_pairs):
        tree_list = sim_coal(num_leaves, 2)  # Simulate a pair of trees instead
        rank_count = bfs_path_rank_sequence(tree_list.trees[0],
                                            tree_list.trees[1])
        print(rank_count)


# use own implementation of coalescent to plot ranked SPR distances (HSPR if hspr=0, otherwise (default) RSPR) between coalescent trees (i.e. uniform ranked trees)
def coal_pw_spr_dist(num_leaves,
                     num_tree_pairs,
                     hspr=False,
                     output_file='',
                     distances_file=''):
    # Plotting the distances for num_tree_pairs simulated pairs of trees and save plot (if filehandle given) in output_file
    distances = []

    for i in range(0, int(num_tree_pairs)):
        if i % 100 == 0:
            print('iteration', i)
        # Simulate a pair of trees instead of a list with num_tree trees
        tree_list = sim_coal(num_leaves, 2)
        distances.append(
            len(rankedspr_bfs(tree_list.trees[0], tree_list.trees[1])) - 1)
    if distances_file != '':
        np.savetxt(distances_file, distances, delimiter=' ')
    # Plot histogram
    d = pd.DataFrame(data=distances)
    upper_bound = max(distances)
    b = np.arange(-.5, upper_bound + 1.5, 1)
    sns.set_theme(font_scale=1.2)
    sns.histplot(d,
                 color='#b02538',
                 edgecolor='black',
                 alpha=1,
                 binwidth=1,
                 binrange=[-.5, upper_bound + 1.5],
                 stat='density',
                 legend=False)
    plt.xlabel("Distance")
    plt.ylabel("Proportion of trees")
    if hspr == False:
        plt.savefig("output/plots/rspr_distribution_" + str(num_leaves) +
                    "_n_" + str(num_tree_pairs) + ".eps")
    else:
        plt.savefig("output/plots/hspr_distribution_" + str(num_leaves) +
                    "_n_" + str(num_tree_pairs) + ".eps")
    plt.clf()
    # plt.show()
    # plts.plot_hist(distances, bins, output_file)


# use own implementation of coalescent to plot difference in RSPR/HSPR distances between two coalescent trees and the same trees with one leaf deleted
def distance_del_leaf(num_leaves,
                      num_deletions,
                      num_tree_pairs,
                      hspr=False,
                      output_file='',
                      distances_file=''):
    # Plotting the distances for num_tree_pairs simulated pairs of trees and save plot (if filehandle given) in output_file
    distances = []
    if num_leaves <= 6:
        (dist, tree_dict,
         tree_index_dict) = get_distance_matrix(num_leaves, hspr)

    for i in range(0, int(num_tree_pairs)):
        # Simulate a pair of trees instead of a list with num_tree trees
        tree_list = sim_coal(num_leaves, 2)
        if num_leaves <= 6:
            tree_str1 = tree_to_cluster_string(tree_list.trees[0])
            tree_str2 = tree_to_cluster_string(tree_list.trees[1])
            d = dist[tree_dict[tree_str1]][tree_dict[tree_str2]]
        else:
            d = len(rankedspr_bfs(tree_list.trees[0], tree_list.trees[1])) - 1

        # # try deleting every leaf and see how distance decreases
        # max_dist = d
        # for i in range(0,num_leaves):
        #     tree1 = del_leaf(tree_list.trees[0],i)
        #     tree2 = del_leaf(tree_list.trees[1],i)
        #     current_dist = len(rankedspr_bfs(tree1, tree2))-1
        #     if current_dist<max_dist:
        #         max_dist = current_dist
        # distances.append(d-max_dist)

        # # alternatively: try to delete every pair of leaves and look at minimum distance
        # current_dist = []
        # for i in range(0,num_leaves-1):
        #     tree1 = del_leaf(tree_list.trees[0],i)
        #     tree2 = del_leaf(tree_list.trees[1],i)
        #     for j in range(0,num_leaves-2):
        #         tree1_1 = del_leaf(tree1,j)
        #         tree2_1 = del_leaf(tree2,j)
        #         current_dist.append(len(rankedspr_bfs(tree1_1, tree2_1))-1)
        #     print(d - max(current_dist), i, j, tree_to_cluster_string(tree_list.trees[0]), tree_to_cluster_string(tree_list.trees[1]))
        # distances.append(d - max(current_dist))

        # even another alternative: delete the two cherry leaves
        tree1 = tree_list.trees[0]
        tree2 = tree_list.trees[1]
        c1 = min(tree1.node_array[num_leaves].children[0],
                 tree1.node_array[num_leaves].children[1])
        c2 = max(tree1.node_array[num_leaves].children[0],
                 tree1.node_array[num_leaves].children[1])
        tree1 = del_leaf(tree1, c2)
        tree1 = del_leaf(tree1, c1)
        tree2 = del_leaf(tree2, c2)
        tree2 = del_leaf(tree2, c1)
        d1 = len(rankedspr_bfs(tree1, tree2)) - 1
        distances.append(d - d1)
        if d == math.floor(3 / 2 * (num_leaves - 2)):
            print("Diameter distance", math.floor(3 / 2 * (num_leaves - 2)),
                  " Difference in distances after deleting cherry:", d - d1)

        # if d-d1 == 3:
        #     print("original trees:")
        #     print(tree_to_cluster_string(tree_list.trees[0]))
        #     print(tree_to_cluster_string(tree_list.trees[1]))
        #     print("trees after deleting leaves:")
        #     print(tree_to_cluster_string(tree1))
        #     print(tree_to_cluster_string(tree2))
    print(distances)

    print("maximum differences in distances:", max(distances))
    if distances_file != '':
        np.savetxt(distances_file, distances, delimiter=' ')
    # Plot histogram
    d = pd.DataFrame(data=distances)
    upper_bound = max(distances)
    b = np.arange(-.5, upper_bound + 1.5, 1)
    sns.set_theme(font_scale=1.2)
    sns.histplot(d,
                 color='#b02538',
                 edgecolor='black',
                 alpha=1,
                 binwidth=1,
                 binrange=[-.5, upper_bound + 1.5],
                 stat='density',
                 legend=False)
    plt.xlabel("Distance")
    plt.ylabel("Proportion of trees")
    if hspr == False:
        plt.savefig("output/plots/rspr_dist_diff_" + str(num_leaves) + "_n_" +
                    str(num_tree_pairs) + ".eps")
    else:
        plt.savefig("output/plots/hspr_dist_diff_" + str(num_leaves) + "_n_" +
                    str(num_tree_pairs) + ".eps")
    plt.clf()
    # plt.show()
    # plts.plot_hist(distances, bins, output_file)


# check for all pairs of trees on num_leaves leaves how the distance changes when deleting a leaf (for evey leaf)
# only works for small number of leaves, because we need the distance matrix!
def full_distance_del_leaf(num_leaves, hspr=False, distances_file=''):
    # Plotting the distances for num_tree_pairs simulated pairs of trees and save plot (if filehandle given) in output_file
    distances = []
    (dist, tree_dict, tree_index_dict) = get_distance_matrix(num_leaves, hspr)
    (dist_small, tree_dict_small,
     tree_index_dict_small) = get_distance_matrix(num_leaves - 1, hspr)
    tree_pairs_increasing_distance = dict()
    increasing_dist_tp = dict()
    num_increasing_tp = 0

    for i in range(0, len(dist)):
        tree1_str = tree_index_dict[i]
        tree1 = read_from_cluster(tree1_str)
        for j in range(i + 1, len(dist)):
            tree2_str = tree_index_dict[j]
            tree2 = read_from_cluster(tree2_str)
            # because of symmerty is sufficient to consider only cases where leaf at fixed position (num_leaves-1) is deleted
            tree1_d = del_leaf(tree1, num_leaves - 1)
            tree2_d = del_leaf(tree2, num_leaves - 1)
            tree1_d_str = tree_to_cluster_string(tree1_d)
            tree2_d_str = tree_to_cluster_string(tree2_d)
            small_distance = dist_small[tree_dict_small[tree1_d_str]][
                tree_dict_small[tree2_d_str]]
            distances.append(dist[i][j] - small_distance)
            if dist[i][j] - small_distance < 0:
                if small_distance - dist[i][j] not in increasing_dist_tp:
                    increasing_dist_tp[small_distance - dist[i][j]] = []
                increasing_dist_tp[small_distance - dist[i]
                                   [j]].append(set([tree1_str, tree2_str]))
                num_increasing_tp += 1
                if small_distance - dist[i][
                        j] not in tree_pairs_increasing_distance:
                    tree_pairs_increasing_distance[small_distance -
                                                   dist[i][j]] = []
                    if small_distance != rnni_distance(tree1_d, tree2_d):
                        print("small distance != RNNI distance for")
                        print(tree1_d_str, tree2_d_str)
                add_unique_topology_pair(
                    tree_pairs_increasing_distance[small_distance - dist[i][j]], tree1, tree2)

    # Print which tree pairs decrease distance and by how much
    for d in tree_pairs_increasing_distance:
        print("Tree pairs with distance", d)
        for tp in tree_pairs_increasing_distance[d]:
            print(tp)
    for i in increasing_dist_tp.keys():
        print("For", len(increasing_dist_tp[i]),
              "tree pairs the distance increases by", i,
              "when deleting one leaf.")
    if distances_file != '':
        np.savetxt(distances_file, distances, delimiter=' ')

    # Plot histogram
    if hspr == True:
        f = "output/plots/full_hspr_dist_diff_" + str(num_leaves) + "_n.eps"
    else:
        f = "output/plots/full_rspr_dist_diff_" + str(num_leaves) + "_n.eps"

    plot_hist(distances,
              xlabel='difference in distance',
              ylabel='count',
              filehandle=f,
              density=False)


def rank_moves_distribution(num_leaves):
    (d, tree_dict, tree_index_dict) = get_distance_matrix(num_leaves,
                                                           hspr=False)
    max_dist = np.amax(d)
    rank_move_dict = dict()  # keys: distances between trees, values: lists of numbers of rank moves on every shortest path between every pair of trees with corresponding distance
    # initialise rank_move_dict:
    for i in range(1, max_dist + 1):
        rank_move_dict["dist" + str(i)] = [0] * max_dist

    # for every pair of trees, add list with number of rank moves on shortest paths to rank_move_dict[d] where d is the distance between them
    for i in range(0, len(d)):
        # if i % (math.floor(len(d) / 100)) == 0:
        #     print("progress:", int(100 * i / len(d)), "percent")
        tree1_str = tree_index_dict[i]
        tree1 = read_from_cluster(tree1_str)
        for j in range(i + 1, len(d)):
            tree2_str = tree_index_dict[j]
            tree2 = read_from_cluster(tree2_str)
            # list of number of rank moves on all shortest paths between tree1 and tree2
            rm = count_rank_moves_all_shortest_paths(
                tree1, tree2, d, tree_dict, tree_index_dict)
            for k in range(0, max(rm) + 1):
                rank_move_dict["dist" + str(d[i][j])][k] += rm.count(k)

    plt.clf()
    # Plot number of rank moves per shortest paths, one line for each possible distance
    d = pd.DataFrame(data=rank_move_dict)
    print(d)
    sns.set_theme(font_scale=2, style='whitegrid')
    sns.lineplot(data=d, markers=True)
    plt.xlabel("Number of rank moves on shortest path")
    plt.ylabel("Number of paths")
    plt.savefig("output/plots/rank_move_distribution_" +
                str(num_leaves) + "_n.eps")
    # plt.clf()
    # plt.show()

    # also plot the 'normalised' number of rank moves on shortest path, i.e. divide the number of paths with x rank moves by the total number of paths:
    norm = dict()
    for i in rank_move_dict:
        norm[i] = [float(x / sum(rank_move_dict[i]))
                   for x in rank_move_dict[i]]
    print(norm)

    # Plot relative number of rank moves per shortest paths, one line for each possible distance
    plt.clf()
    norm_d = pd.DataFrame(data=norm)
    print(norm_d)
    sns.set_theme(font_scale=2, style="whitegrid")
    sns.boxplot(data=norm_d)
    sns.stripplot(data=norm_d, size=4)
    plt.xlabel("Length of shortest paths")
    plt.ylabel("Relative number of rank moves")
    plt.tight_layout()
    plt.savefig("output/plots/rank_move_distribution_norm_" + str(num_leaves) +
                "_n_boxplot.eps")
    # plt.show()

    # also do a lineplot for normalised number of rank moves on paths
    plt.clf()
    d = pd.DataFrame(data=norm)
    sns.set_theme(font_scale=2, style="white")
    sns.lineplot(data=d, markers=True)
    # plt.xlabel("Normalised number of rank moves on shortest path")
    # plt.ylabel("Number of paths")
    # # plt.tight_layout()
    plt.savefig("output/plots/rank_move_distribution_norm_" + str(num_leaves) +
                "_n.eps")


# Compute all shortest path in HSPR (using the distance matrix for the whole tree space computed by SEIDEL)
# Does not actually output those shortest paths, but prints predecessor list
def all_shortest_paths(tree1, tree2):
    num_leaves = tree1.num_leaves
    if not exists("output/distance_matrix_" + str(num_leaves) + "_leaves_hspr.npy"):
        rankedspr_seidel(num_leaves, hspr = True)
    (d, tree_dict, tree_index_dict) = get_distance_matrix(num_leaves, True)
    tree1_str = tree_to_cluster_string(tree1)
    tree2_str = tree_to_cluster_string(tree2)

    tree1_index = tree_dict[tree1_str]
    tree2_index = tree_dict[tree2_str]

    distance = d[tree1_index][tree2_index]

    # for every tree that is on a shortest path, save all predecessors of it in dictionary pred:
    pred = dict()
    for tree_index in range(0, len(d)):
        if d[tree1_index][tree_index] + d[tree_index][tree2_index] == distance:
            tree = read_from_cluster(tree_index_dict[tree_index])
            neighbourhood = hspr_neighbourhood(tree)
            for i in range(0, neighbourhood.num_trees):
                predecessor = neighbourhood.trees[i]
                pred_str = tree_to_cluster_string(predecessor)
                pred_index = tree_dict[pred_str]
                if d[tree1_index][pred_index] + d[pred_index][tree_index] + d[
                        tree_index][
                            tree2_index] == distance:  # if predecessor is on shortest path from tree1 to tree2
                    if tree_index in pred:
                        pred[tree_index].add(pred_index)
                    else:
                        pred[tree_index] = set([pred_index])

    for i in pred:
        print("predecessors of ",tree_index_dict[i])
        for k in pred[i]:
            print(tree_index_dict[k])
