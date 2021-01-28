# RNNI distance Distributions

## Simulate trees

### Coalescent
- Simulating coalescent trees using R package ape (rcoal) with script *src/simulate_coalescent.R*
- Saved to file *simulations/simulated_trees/coal/y/coal_trees_x_n.nex* where x number of leaves and y number of trees

### Birth-Death
- Simulating birth-death trees using R package geiger (sim.bdtree) with script *src/simulate_bd.R*
- saved in *simulations/simulated_trees/bd/y/bd_trees_x_n.nex* where x number of leaves and y number of trees


## Analysing distance distributions

- Plotting distances between pairs of trees with index i,i+1 for even i (from trees simulated as above), using *src/analyse_distance_distributions.py* (rnni_distances_tree_pairs/ rf_distances_tree_pairs)
- Plots are saved as *simulations/distance_distributions/D/rnni_distribution_x_n_y_m.eps* (for RNNI) and *simulations/distance_distributions/D/rf_distribution_x_n_y_N.eps* (for RF), where D is the prior distribution (coalescent or bd)

### Distances from focal tree

- Plotting distances between tree that is chosen randomly at uniform to all other trees (from trees simulated as above), using *src/analyse_distance_distributions.py*
- Plots are saved as *simulations/distance_distributions/D/rnni_focal_distribution_x_n_y_m.eps* (for RNNI) and *simulations/distance_distributions/D/rf_focal_distribution_x_n_y_N.eps* (for RF), where D is the prior distribution (coalescent or bd)

### Mean distance vs. n

- Plotting mean distance (mean of the ones produced by rnni_distances_tree_pairs as described above) vs. number of leaves (n=1:40) -- saved as *simulations/distance_distribution/D/rnni_meanVSn_x_n_y_N*, where D is the prior distribution (coalescent or bd)


## Simulate posterior sample

### Alignments

(1) Simulate **generating trees** with R: rcoal/sim.bdtree (coalescent/birth-death) and save in *simulations/posterior/bd/D/*, where D is the prior (coal or bd). Name of files: *coal_20_leaves.new* / *bd_tree_b1_d0_20_leaves.new* (coal/bd)  
Resulting trees have 20 leaves and alignments have length 10,000
(2) Simulate alignments along those trees using simSeq and JC model. Alignments are saved in same folder as trees as *coal_alignment_20_sequences_10000_length.nex* / *bd_alignment_b1_d0_20_sequences_10000_length.nex* (coal/bd).
(3) Run Beast using the xml file (produced by Beauti) in the same folder as generating tree and alignment  
Settings for beast: JC model, strict clock, coal/bd prior, 1,000,000 trees, saving every 1,000st tree.
(4) Analyse pw distances between each tree and the following one (d(T_i,T_{i+1}) for all i) with function *rnni_distances_consecutive_tree_pairs* and plot distances to *rnni_consecutive_pairs_dist.eps* in the same folder