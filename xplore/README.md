# RNNI distance Distributions

## Simulate trees

### Coalescent
- Simulating coalescent trees using R package ape (rcoal) with script "./src/simulate_coalescent.R"
- Saved to file "./simulated_trees/coal/y/coal_trees_x_n.nex" where x number of leaves and y number of trees

### Birth-Death
- Simulating birth-death trees using R package geiger (sim.bdtree) with script "./src/simulate_bd.R"
- saved in "simulated_trees/bd/y/bd_trees_x_n.nex" where x number of leaves and y number of trees


## Analysing distribution

- Plotting distances between pairs of trees with index i,i+1 for even i (from trees simulated as above), using "./src/analyse_distance_distributions.py" (rnni_distances_tree_pairs/ rf_distances_tree_pairs)
- Plots are saved as "./distance_distributions/D/rnni_distribution_x_n_y_m.eps" (for RNNI) and "./distance_distributions/D/rf_distribution_x_n_y_N.eps" (for RF), where D is the prior distribution (coalescent or bd)

### Distances from focal tree

- Plotting distances between tree that is chosen radnomly at uniform to all other trees (from trees simulated as above), using "./src/analyse_distance_distributions.py"
- Plots are saved as "./distance_distributions/D/rnni_focal_distribution_x_n_y_m.eps" (for RNNI) and "./distance_distributions/D/rf_focal_distribution_x_n_y_N.eps" (for RF), where D is the prior distribution (coalescent or bd)

### Mean distance vs. n

- Plotting mean distance (mean of the ones produced by rnni_distances_tree_pairs as described above) vs. number of leaves (n=1:40) -- saved as "./distance_distribution/D/rnni_meanVSn_x_n_y_N", where D is the prior distribution (coalescent or bd)