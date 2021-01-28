# Author: Lena Collienne
# simulate coalescent trees

library(geiger)
library(rstudioapi)

# Setting the working directory -- might need to be changes when not using RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../simulated_trees")

bd_tree_sim <- function(num_leaves, num_trees) {
  trees = c()
  for (i in 1:num_trees){
    trees[[i]] <- sim.bdtree(b=1, d=0, stop="taxa", n=num_leaves)
  }
  write.nexus(trees, file = paste("bd/", num_trees, "/bd_trees_", num_leaves, "_n.nex", sep = ""))
}

# bd_tree_sim(100,20000)

for(num_leaves in 3:100){
  bd_tree_sim(num_leaves,2000)
}