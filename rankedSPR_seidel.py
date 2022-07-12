__author__ = 'Lena Collienne'
# Compute distance matrix for rankedSPR graph (from adjacency matrix, using SEIDEL implementation from RNNI_code package)
import sys
sys.path.append('treeOclock/')
sys.path.append('treeOclock/dct_parser/')

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import time
from rankedSPR import *
from os.path import exists


_seidel = ctypes.CDLL("./libseidel.so")
_seidel.test_function.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32)
_seidel.seidel.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32)
_seidel.seidel_recursive.argtypes = (ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int32, ctypes.c_int32)

def rankedspr_seidel(n, hspr=1):
    # compute distance matrix for RSPR (or HSPR if HSPR=0), for trees on n leaves
    print('number of leaves:', n)
    AI = rankedSPR_adjacency(n, hspr)
    A = np.ascontiguousarray(AI[0], dtype=np.int32)
    time1 = time.time()
    _seidel.seidel(A, A.shape[0])
    if(hspr==0):
        np.save('SPR/distance_matrix_' + str(n) + '_leaves_hspr', A)
    else:
        np.save('SPR/distance_matrix_' + str(n) + '_leaves', A)
    time2 = time.time()
    print("C Seidel took {:.3f}ms".format((time2 - time1)*1000.0))
    print("diameter: ", np.amax(A))

def rankedspr_wo_RNNI_seidel(n):
    # compute distance matrix for RSPR (or HSPR if HSPR=0), for trees on n leaves
    print('number of leaves:', n)
    AI = rankedSPR_wo_RNNI_adjacency(n)
    A = np.ascontiguousarray(AI[0], dtype=np.int32)
    time1 = time.time()
    _seidel.seidel(A, A.shape[0])
    np.save('SPR/wo_RNNI_distance_matrix_' + str(n) + '_leaves', A)
    time2 = time.time()
    print("C Seidel took {:.3f}ms".format((time2 - time1)*1000.0))
    print("diameter: ", np.amax(A))

