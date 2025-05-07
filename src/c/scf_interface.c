#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// Fortran subroutine interface
void run_scf_interface_py(int32_t *natoms, int32_t *nelectrons, int32_t *atomic_numbers, double *coords, 
                          double *eps_scf, int32_t *max_iter, double *final_energy, 
                          int32_t *nbasis_atom, double *basis_exp, int32_t *basis_set);

void run_scf_c(int natoms, int nelectrons, int *atomic_numbers, double *coords, double eps_scf, 
               int max_iter, double *final_energy, int *nbasis_atom, double *basis_exp, int basis_set) {
    // Call the Fortran subroutine
    run_scf_interface_py(&natoms, &nelectrons, atomic_numbers, coords, &eps_scf, &max_iter, 
                         final_energy, nbasis_atom, basis_exp, &basis_set);
}
