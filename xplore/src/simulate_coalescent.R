# Author: Lena Collienne
# simulate coalescent trees

library(ape)
library(rstudioapi)

# Setting the working directory -- might need to be changes when not using RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../simulated_trees")

coalescent_tree_sim <- function(num_leaves, num_trees) {
  trees = c()
  for (i in 1:num_trees){
    trees[[i]] <- rcoal(num_leaves)
  }
  write.nexus(trees, file = paste("coal/coal_trees_", num_leaves, "_n_", num_trees, "_N.nex", sep = ""))
}

# coalescent_tree_sim(50,20000)

for(num_leaves in 3:40){
  coalescent_tree_sim(num_leaves,20000)
}