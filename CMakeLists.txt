cmake_minimum_required(VERSION 3.15)
project(scf_project LANGUAGES C Fortran)

# Set the output directory for libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Enable Fortran preprocessing if needed
set(CMAKE_Fortran_PREPROCESS ON)

# Create the Fortran library
add_library(scf_lib
    src/fortran/molecule_module.f90
    src/fortran/integral_module.f90
    src/fortran/utils_module.f90
    src/fortran/scf_module.f90
    src/fortran/scf_interface.f90
    src/fortran/slater.f90
)

# Create the C interface
add_library(scf SHARED src/c/scf_interface.c)

# Link Fortran and C libraries
target_link_libraries(scf PRIVATE scf_lib)

# Link to LAPACK if needed
find_package(LAPACK REQUIRED)
target_link_libraries(scf PRIVATE ${LAPACK_LIBRARIES})

# Install
install(TARGETS scf LIBRARY DESTINATION lib)
