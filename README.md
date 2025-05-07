# CPU/GPU acceleration of SCF module
This repository contains the code and theoretical background of the Project in Practice with Molecular Quantum Solutions produced during the Spring block of 2025.

> Contact:
> Miguel Alonso (Owner): mi.alonso.mediavilla@gmail.com

---

## Input files
The input files are written in the following way:
```console
2 4 1 6 -> n_atom n_electron spin_multiplicity n_total_basis_sto
0.0 0.0 0.0 1.0 2 -> x y z atomic_num n_basis_sto
1.3 -> exponent_sto_1
0.7 -> exponent_sto_2
0.0 0.0 3.1 3.0 4 -> x y z atomic_num n_basis_sto
3.5 -> .
2.0 -> .
0.7 -> .
0.3
```
This is an example of an input file for LiH. 


## Build and usage
The Fortran and C modules can be built by typing in the terminal:
```console
~$ mkdir build
cd build
cmake ..
cmake --build .
```
From there on, one can use the test_scf.ipynb Jupiter notebook to run SCF simulations on the input files.
