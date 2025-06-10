# CPU/GPU acceleration of SCF module
This repository contains the code and theoretical background of the Project in Practice with Molecular Quantum Solutions (MQS) produced during the Spring block of 2025.

---
> Contact:
> Miguel Alonso (Owner): mi.alonso.mediavilla@gmail.com

---

## Features
- Usage of STO-{2-6}g basis sets.
- .XYZ input file support.
- Calculation of SCF energy for unrestricted (open-shell) and restricted (closed-shell) systems.
- Calculation of two-Electron Repulsive Integrals (ERIs) with Libint.
- OpenMP accelerated ERI matrix building.
- OpenACC accelerated Fock matrix building and SCF routine.


## External Dependencies
- Libint.
- Boost (Libint dependency).


## Libint quick build
For an in-depth explanation, visit the [LIBINT repository](https://github.com/evaleev/libint). For convenience, the following can be used with this project:
Download the latest version from (https://github.com/evaleev/libint/releases/tag/v2.11.1).
```console
~$ tar -xvzf libint-2.x.y.tgz
cd libint-2.x.y
cmake . -DCMAKE_INSTALL_PREFIX=../libint-2.11.1-install/ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="-O1" -DREQUIRE_CXX_API=OFF -DENABLE_FORTRAN=ON -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.1/compilers/bin/nvfortran -DCMAKE_POSITION_INDEPENDENT_CODE=ON
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
- Usage of cuSOLVER for eigenvalue-like calculations.
- Inclusion of Libint one-electron integral calculations (currently own implementation).
- Apply OpenMP to the formation of the one-electron integral matrices.
- Benchmarking.
