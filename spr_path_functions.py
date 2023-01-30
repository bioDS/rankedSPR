__author__ = 'Lena Collienne'

import sys
from ctypes import *
sys.path.append('treeOclock')
from tree_functions import *

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/spr_path.so')


rankedspr_path_bottom_up_hspr = lib.rankedspr_path_bottom_up_hspr
rankedspr_path_bottom_up_hspr.argtypes = [POINTER(TREE), POINTER(TREE)]
rankedspr_path_bottom_up_hspr.restype = TREE_ARRAY
