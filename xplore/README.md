# RNNI distance Distributions


## Simulate trees

### Coalescent
- Simulating coalescent trees using R package ape (rcoal) with script *src/simulate_coalescent.R*
- Saved to file *simulations/simulated_trees/coal/y/coal_trees_x_n.nex* where x number of leaves and y number of trees

### Birth-Death
- Simulating birth-death trees using R package geiger (sim.bdtree) with script *src/simulate_bd.R*
- saved in *simulations/simulated_trees/bd/y/bd_trees_x_n.nex* where x number of leaves and y number of trees

### simulate_trees.py
*src/simulate_trees.py* contains simulations for ranked trees as C data structure

  | function			|	Description
  ---    |   ---
  | *sim_coal* | coalescent simulation in |
  | *sim_cat* | simulate caterpillar tree uniformly among caterpillar trees |


## Analysing distance distributions
*src/analyse_distance_distribution.py*

### Using coalescent implementation sim_coal
| Function | Description | Save File for Plot (if existing)
--- | --- | ---
| *coal_pw_dist* | plot distances between tree pairs i,i+1 for even i (distances are independent) in list of trees simulated using *sim_coal* | *simulations/distance_distribution/coalescent/own_coal_distr_20_n_20000_N.eps* |
| *given_focal_distribution* | plot distances from num_tree trees (sim_coal) to given focal tree | for caterpillar tree (sim.sim_cat): *simulations/distance_distribution/coalescent/caterpillar_distances_20_n_20000_N.eps*, for fully balanced tree: *simulations/distance_distribution/coalescent/dist_distribution_to_fully_balanced_16_n_10000_N.eps*, for coalescent tree (uniform): *simulations/distance_distribution/coalescent/coal_focal_dist_20_n_20000_N.eps* |
| *dist_distribution_to_caterpillar* | plot distances between num_tree caterpillar trees (sim_cat) to a fixed coalescent tree | *simulations/distance_distribution/coalescent/dist_distribution_to_caterpillars_20_n_10000_N.eps* |
| *dist_distribution_btw_caterpillar* | plot distances between num_tree caterpillar tree pairs i,i+1 for even i (distances are independent) in list of trees simulated using *sim_cat* | *simulations/distance_distribution/coalescent/dist_distribution_btw_caterpillars_20_n_10000_N.eps* |
| *compare_given_focal_tree_dist* | compare distances between num_trees coalescent trees to fully balanced tree with distances to caterpillar tree | *../simulations/distance_distribution/coalescent/compare_cat_balanced_16_n_10000_N.eps* |
| *compare_expected_dist_to_simulation* | Compare approximation of expected distance between two coalescent trees with mean of simulated distances | *../simulations/distance_distribution/coalescent/compare_expected_dist_to_simulation_40_n_20000_N.eps* |
| *plot_approx_exp_dist* | Plot approximated distance depending on n (input is max n, computes it for n=3:max_num_leaves) | *../simulations/distance_distribution/coalescent/approx_exp_dist_n_3_to_500.eps* |
| *random_walk_distance* | Plots distance between simulated start tree and tree after random walk of length k, num_iterations repetitions of this | *../simulations/distance_distribution/coalescent/random_walk_mean_dist_n_6_k_1_to_1000_N_1000.eps* |
  |  |  
| *mean_distances_n* | takes as input one of the functions above and compute distances according to them for n between given boundaries and plots mean and variance of distances depending on n |*simulations/distance_distribution/coalescent/mean_and_var_dist_n_3_to_40_N_20000.eps* (alternatively also *cat_mean_and_var_dist_n_3_to_40_N_20000.eps*, *btw_cat_mean_and_var_dist_n_3_to_40_N_20000.eps*, *from_cat_mean_and_var_dist_n_3_to_40_N_20000.eps* etc.) |
| *mean_distances_repeat* | takes as input one of the functions above and compute distances according to them, repeating it for a given number of times (num_iterations) and plots mean and variance of distances for each iteration |*simulations/distance_distribution/coalescent/mean_and_var_dist_n_3_to_40_N_20000.eps* (alternatively also *to_cat_mean_and_var_dist_n_3_to_40_N_20000.eps* etc.) |
| *mean_distances_log_n* | takes as input one of the functions above and compute distances according to them for n between 2^i for i=2^2,..2^max_exp (max_exp is input) and plots mean of distances depending on n |*simulations/distance_distribution/coalescent/mean_dist_log_n_10_N_10000.eps* |
| *mean_distances_log_n* | takes as input one of the functions above and compute distances according to them for n between 100*min_num_leaves and 100*max_num_leaves (min_num_leaves and mex_num_leaves are input) and plots mean of distances depending on n |*simulations/distance_distribution/coalescent/mean_dist_100_n_4_to_20_N_1000.eps* |
| *random_walk_mean_distance* | Plots distance between simulated start tree and tree after random walk of length k for varying k, num_iterations repetitions for each k, fixed n; For median=TRUE, return approximation of value the distance between start and end tree converges to | *../simulations/distance_distribution/coalescent/random_walk_dist_n_6_k_20_N_1000.eps* |

### Distance Distributions

- Plotting distances between pairs of trees with index i,i+1 for even i (from trees simulated as above), using *src/analyse_distance_distributions.py* (rnni_distances_tree_pairs/ rf_distances_tree_pairs)
- Plots are saved as *simulations/distance_distributions/D/rnni_distribution_x_n_y_m.eps* (for RNNI) and *simulations/distance_distributions/D/rf_distribution_x_n_y_N.eps* (for RF), where D is the prior distribution (coalescent or bd)
- For own coalescent simulation (*sim_coal*): 

### Distances from random focal tree

- Plotting distances between tree that is chosen randomly at uniform to all other trees (from trees simulated as above), using *src/analyse_distance_distributions.py*
- Plots are saved as *simulations/distance_distributions/D/rnni_focal_distribution_x_n_y_m.eps* (for RNNI) and *simulations/distance_distributions/D/rf_focal_distribution_x_n_y_N.eps* (for RF), where D is the prior distribution (coalescent or bd)

### Mean distance vs. n

- Plotting mean distance (mean of the ones produced by rnni_distances_tree_pairs as described above) vs. number of leaves (n=1:40) -- saved as *simulations/distance_distribution/D/rnni_meanVSn_x_n_y_N*, where D is the prior distribution (coalescent or bd)

### Uniform distances

- exact distance distributions on up to 6 leaves, using the RNNI_code repo (test_distance_distribution)
- distance distributions saved in *distance_distribution/uniform/all_rnni_dist_x_n.eps*, where x is the number of leaves (x=4,5,6)


## Posterior sample

### Simulation

(1) Simulate **generating trees** with R: rcoal/sim.bdtree (coalescent/birth-death) and save in *simulations/posterior/bd/D/*, where D is the prior (coal or bd). Name of files: *coal_20_leaves.new* / *bd_tree_b1_d0_20_leaves.new* (coal/bd)
Resulting trees have 20 leaves and alignments have length 10,000  
(2) Simulate alignments along those trees using simSeq and JC model. Alignments are saved in same folder as trees as *coal_alignment_20_sequences_10000_length.nex* / *bd_alignment_b1_d0_20_sequences_10000_length.nex* (coal/bd).  
(3) Run Beast using the xml file (produced by Beauti) in the same folder as generating tree and alignment  
Settings for beast: JC model, strict clock, coal/bd prior, 1,000,000 trees, saving every 1,000st tree.

### Analysis

- We do currently not discard a burn-in!
- Analyse pw distances between each tree and the following one (d(T_i,T_{i+1}) for all i) with function *rnni_distances_consecutive_tree_pairs* and plot distances to *rnni_consecutive_pairs_dist.eps* in the same folder  
This can also been done with the RF distance (functions in *src/rf_distance_analysis.py*)  
- Plot hist of all pw distances (RNNI and RF) with pw_rnni_dist/pw_rf_dist -- output in (*simulations/posterior/D* for distribution D=coal/bd) as *rnni_all_pw_dist.eps*/*rf_all_pw_dist.eps*  
  This is to get an impression of how 'dense' the chain is for the corresponding distance measure (RNNI vs RF)
- Plot all distances from trees in posterior sample to mcc tree (*simulations/posterior/D/mcc_summary.tree*, computed with treeannotator) as *../simulations/posterior/D/z_mcc_dist.eps*, where z is rnni or rf, and D the distribution bd or coal.