# PyFAcc-SCF
This repository contains the Python Fortran Accelerated Self-Consistent Field (PyFAcc-SCF) project. This repository was created as part of a Project in Practice and a later MSc Thesis project at the University of Copenhagen in collaboration with Molecular Quantum Solutions (MQS) from April 2025 to May 2026.


---
> Contact:
> Miguel Alonso (Owner): mi.alonso.mediavilla@gmail.com

---
## Objectives
Creation of a multi-language GPU-accelerated SCF routine for the study of molecular systems. Usage of OpenACC/OpenMP/OpenMPI to adapt to heterogeneous systems and different accelerators. 

## Current Features
- Usage of STO-{2-6}g basis sets.
- .XYZ input file support.
- Calculation of SCF energy for unrestricted (open-shell) and restricted (closed-shell) systems.
- Calculation of two-electron Repulsive Integrals (ERIs) with Libint.
- OpenACC accelerated one-electron integrals, Fock matrix building, and SCF energy calculation.


## External Dependencies
- Libint.
- Boost (Libint dependency).
- Eigen3 (Libint dependency).


## Libint quick build
For an in-depth explanation, visit the [LIBINT repository](https://github.com/evaleev/libint). 
For convenience, the following can be used with this project:

Download the latest version from (https://github.com/evaleev/libint/releases/tag/v2.11.1).
```console
~$ tar -xvzf libint-2.x.y.tgz
cd libint-2.x.y
cmake . -DCMAKE_INSTALL_PREFIX=../libint-2.11.1-install/ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="-O1" -DENABLE_FORTRAN=ON -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.1/compilers/bin/nvfortran -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build .
cmake --build . --target install
```
Make sure to change the path of installation of the compiler (nvfortran) and libint-2.x.y to the downloaded version.



## Build and Usage
The Fortran and C modules can be built by:
```console
~$ mkdir build
cd build
cmake ..
cmake --build .
```
From there on, one can use the included Python file to run SCF simulations on the input files.

---

## Future Improvements and Next Steps
- Addition of Schwarz screening to ERI computations.
- Focus on ERI computations: research on advanced techniques, higher GPU offloading. HGP algorithm for ERI computation following ["Advanced Techniques for High-Performance Fock Matrix Construction on GPU Clusters"](https://pubs.acs.org/doi/10.1021/acs.jctc.4c00994).
- Multi-GPU support via MPI.
- Possible usage of cuSOLVER for eigenvalue-like and linear algebra calculations (fully CUDA package from Nvidia).
- Expansion to other Basis Sets (e.g., cc-pVDZ, cc-pVTZ).
