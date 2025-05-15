cmake_minimum_required(VERSION 3.15)
project(scf_project LANGUAGES C Fortran)

# Set compilers
set(CMAKE_Fortran_COMPILER nvfortran)
set(CMAKE_C_COMPILER nvc)

# Enable Fortran preprocessing and position-independent code
set(CMAKE_Fortran_PREPROCESS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set output directory for libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# OpenACC and optimization flags
set(OPENACC_FLAGS -acc -gpu=cuda12.6 -Minfo=acc -V)

# Find LAPACK
find_package(LAPACK REQUIRED)

# Fortran library with OpenACC
add_library(scf_lib
    src/fortran/molecule_module.f90
    src/fortran/integral_module.f90
    src/fortran/utils_module.f90
    src/fortran/scf_module.f90
    src/fortran/scf_interface.f90
    src/fortran/slater.f90
    src/fortran/structs/basis_function_struct.f90
    src/fortran/structs/constants_struct.f90
    src/fortran/structs/molecule_struct.f90
)
target_compile_options(scf_lib PRIVATE ${OPENACC_FLAGS})
target_link_libraries(scf_lib PRIVATE ${LAPACK_LIBRARIES})

# C interface (with OpenACC support)
add_library(scf SHARED src/c/scf_interface.c)
target_compile_options(scf PRIVATE ${OPENACC_FLAGS})
target_link_options(scf PRIVATE -acc -gpu=cuda12.6)
target_link_libraries(scf PRIVATE scf_lib ${LAPACK_LIBRARIES})

# Install target
install(TARGETS scf LIBRARY DESTINATION lib)
