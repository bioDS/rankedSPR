# README

This repo contains computations in ranked SPR spaces HSPR and RSPR as defined in the paper *Ranked Subtree Prune and Regraft* by Lena Collienne, Chris Whidden, and Alex Gavryushkin (2023).


## Compiling

Download repository, get all submodules, and compile C code (using gcc):

    git clone git@github.com:bioDS/rankedSPR.git
    git submodule update --recursive --init
    cd rankedSPR
    make


## Functions executable from Python

There are two functions in *test_for_paper.py* that are used to show statements in the paper.
Firstly, we can show computationally that the diameter for small number ($n \leq 7$) of leaves follows the formula $\lfloor\frac{3}{2} (n-2)\rfloor$:  
> We compute the HSPR distance for all pairs of ranked trees, using Seidel's algorithm.  
> The maximum of those distances is returned by `get_diameter(num_leaves)`.

Furthermore, we show that the algorithm for computing a path as given in the proof of Theorem 1 does not compute shortest paths:  
> This algorithm is implemented as `rankedspr_path_bottom_up_hspr(tree1, tree2)` in *spr_path.c*.  
> The function `check_approx_alg(num_leaves)` in *test_for_paper.py* compares if the lengths of shortest path given by that function are identical to their exact distance, using the distance matrix computed by Seidel's algorithm.  

At last, we check that the two trees `[{4,5}:1,{1,2}:2,{1,2,3}:3,{1,2,3,4,5}:4]` and `[{1,2}:1,{4,5}:2,{1,2,3}:3,{1,2,3,4,5}:4]` have no shortest path between them that contains their shared cluster `{1,2,3}` (Figure 7).
Again, we use Seidel's algorithms to get all shortest paths (see `all_shortest_paths(tree1, tree2)` in *rankedspr_exploration.py*).

Executing `python test_for_paper.py`, gives the following results:

    The diameter of HSPR space for 4 leaves is 3
    For the trees [{1,2}:1,{1,2,3}:2,{1,2,3,4}:3] and [{3,4}:1,{1,3,4}:2,{1,2,3,4}:3] the approximated HSPR distance is 3
    The exact distance is 2
    The shared clusters {1,2,3} is not present in a tree on a shortes path between [{4,5}:1,{1,2}:2,{1,2,3}:3,{1,2,3,4,5}:4] and [{1,2}:1,{4,5}:2,{1,2,3}:3,{1,2,3,4,5}:4]


When changing the number of leaves *num_leaves* in the main function of *test.py*, results for that number of leaves are printed.  
Note that computations beyond 7 leaves might exceed the memory of a laptop, as the distance matrix for all ranked trees on 7 leaves is needed, which requires a lot of memory.