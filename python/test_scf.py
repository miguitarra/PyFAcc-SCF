from parser import parser_input
import ctypes
import numpy as np
import os

# Read input file
input_file = os.path.join(os.path.dirname(__file__), 'h2.in')
natoms, atomic_numbers, coords, nbasis_atom, basis_exp = parser_input(input_file)
print("Number of atoms:")
print(natoms)
print("-" * 40)

print("Atomic numbers:")
print(atomic_numbers)
print("-" * 40)

print("Coordinates:")
print(coords)
print("-" * 40)

print("Number of basis functions per atom:")
print(nbasis_atom)
print("-" * 40)

print("Basis function exponents:")
print(basis_exp)
print("-" * 40)

eps_scf = 1e-6
max_iter = 100
final_energy = ctypes.c_double()

# Load library
lib_path = os.path.join(os.path.dirname(__file__), '..', 'build', 'libscf.dylib')
lib = ctypes.CDLL(lib_path)

# Define prototype
lib.run_scf_interface_py.argtypes = [
    ctypes.c_int, # Number of atoms
    ctypes.POINTER(ctypes.c_int), # Array of atomic numbers
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags="C_CONTIGUOUS"), # Coordinates
    ctypes.c_double, # SCF convergence threshold
    ctypes.c_int, # Maximum number of SCF iterations
    ctypes.POINTER(ctypes.c_double), # Final energy
    ctypes.POINTER(ctypes.c_int), # Number of basis functions per atom
    ctypes.POINTER(ctypes.c_double) # Basis function exponents
]

# Call SCF
lib.run_scf_interface_py(natoms, atomic_numbers, coords, eps_scf, max_iter, ctypes.byref(final_energy), nbasis_atom, basis_exp )

print(f"Final SCF Energy: {final_energy.value:.8f} Hartree")
