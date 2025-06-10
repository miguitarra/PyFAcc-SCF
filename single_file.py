import sys
from sys import platform
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'src', 'python'))
from helper_functions import run_scf
import numpy as np
import matplotlib.pyplot as plt    
import pandas as pd
import ctypes

###################################################################
# Constants
eps_scf = 1e-10
max_iter = 10
final_energy = ctypes.c_double()


###########################################################################

#Library definition:
if platform == 'darwin':
    lib_path = os.path.join(os.path.dirname(__file__), 'build', 'lib', 'libscf.dylib')
else:
    lib_path = os.path.join(os.path.dirname(__file__), 'build', 'lib', 'libscf.so')
lib = ctypes.CDLL(lib_path)

lib.run_scf_interface_py.argtypes = [
    ctypes.c_int,                      # natoms
    ctypes.c_int,                      # nelectrons
    ctypes.c_int,                      # nalpha
    ctypes.c_int,                      # nbeta
    ctypes.c_char,                     # species
    ctypes.POINTER(ctypes.c_double),   # coords
    ctypes.c_double,                   # eps_scf
    ctypes.c_int,                      # max_iter
    ctypes.POINTER(ctypes.c_double),   # final_energy
    ctypes.c_int,                      # basis_set_n
    ctypes.c_int                       # z_num
]


##############################################################################

# Load library
basis = 'sto-3g'
input_file = 'h2o.xyz'

energy = run_scf(lib, input_file, eps_scf, max_iter, final_energy, basis)
print(energy)