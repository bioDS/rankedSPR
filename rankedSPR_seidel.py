__author__ = 'Lena Collienne'
# Compute distance matrix for rankedSPR graph (from adjacency matrix, using SEIDEL implementation from RNNI_code package)
import sys
sys.path.append('treeOclock/')
sys.path.append('treeOclock/dct_parser/')

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import time
from treeOclock.dct_parser.tree_io import *
from treeOclock import *
from rankedSPR import *
from os.path import exists


_seidel = ctypes.CDLL("./libseidel.so")
_seidel.test_function.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32)
_seidel.seidel.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32)
_seidel.seidel_recursive.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32, ctypes.c_int32)


AI = rankedSPR_adjacency(4)
A = np.ascontiguousarray(AI[0], dtype=np.int32)	
time1 = time.time()
_seidel.seidel(A, A.shape[0])
time2 = time.time()
print("C Seidel took {:.3f}ms".format((time2 - time1)*1000.0))
print(A)