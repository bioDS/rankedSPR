__author__ = 'Lena Collienne'
# Testing
import sys
sys.path.append('..')

# from treeOclock.tree_parser.tree_io import *
# from simulate_trees import *
from rankedSPR_seidel import *
from rankedspr_exploration import *

n = 5
d = get_distance_matrix(n, hspr=True)
print("diameter of HSPR for", n, "leaves:", np.amax(d[0]))