cmake_minimum_required(VERSION 3.15)
project(scf_project LANGUAGES C Fortran CXX)

# Set compilers
#set(CMAKE_Fortran_COMPILER nvfortran)
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_C_COMPILER clang)

# Enable Fortran preprocessing and position-independent code
set(CMAKE_Fortran_PREPROCESS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set output directory for libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# OpenMP flags for gfortran
set(OPENMP_FLAGS -fopenmp)

# Find LAPACK
find_package(LAPACK REQUIRED)

# ---- LIBINT CONFIGURATION ----
set(LIBINT2_ROOT ${CMAKE_SOURCE_DIR}/external/libint/libint-2.11.1-install)
set(LIBINT2_INCLUDE_DIR "${LIBINT2_ROOT}/include")
set(LIBINT2_MOD_DIR "${LIBINT2_ROOT}/include")
set(LIBINT2_LIBRARY "${LIBINT2_ROOT}/lib/libint2.a") # Change to .so if you use shared

# ---- END LIBINT CONFIGURATION ----

# Fortran library
add_library(scf_lib
    src/fortran/structs/basis_function_struct.f90  
    src/fortran/structs/constants_struct.f90
    src/fortran/structs/molecule_struct.f90
    src/fortran/structs/eri_struct.f90
    src/fortran/basis.f90
    src/fortran/nuclear.f90
    src/fortran/electronic.f90
    src/fortran/gaussian.f90
    src/fortran/io.f90
    src/fortran/kinetic.f90
    src/fortran/la.f90
    src/fortran/main.f90
    src/fortran/molecule.f90
    src/fortran/overlap.f90
    src/fortran/rscf.f90
    src/fortran/uscf.f90
    src/fortran/utils.f90
    src/fortran/integrals.f90
    src/fortran/libint_f.f90
    src/fortran/boys.f90
)

target_include_directories(scf_lib PRIVATE
    ${LIBINT2_INCLUDE_DIR}
    ${LIBINT2_MOD_DIR}
    ${CMAKE_SOURCE_DIR}/external/libint/libint-2.11.1/fortran
)

target_link_libraries(scf_lib PRIVATE ${LAPACK_LIBRARIES} ${LIBINT2_LIBRARY} stdc++)
target_compile_options(scf_lib PRIVATE ${OPENMP_FLAGS})
target_link_options(scf_lib PRIVATE )

# C/C++ interface
add_library(scf SHARED src/c/scf_interface.c src/cpp/libint_wrapper.cpp) # Use .cpp if using C++ code!
target_link_libraries(scf PRIVATE scf_lib ${LAPACK_LIBRARIES} ${LIBINT2_LIBRARY} stdc++)
target_include_directories(scf PRIVATE ${LIBINT2_INCLUDE_DIR})

# Only add OpenACC flags if using an OpenACC-capable compiler (not gfortran)
if(NVFORTRAN)
    target_link_options(scf PRIVATE -acc -gpu=cuda12.6)
endif()

install(TARGETS scf LIBRARY DESTINATION lib)