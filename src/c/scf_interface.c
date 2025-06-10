#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// Fortran subroutine interface
void run_scf_interface_py(int32_t *natoms, int32_t *nelectrons, int32_t *nalpha, int32_t *nbeta,
                          const char *species, double *coords, double *eps_scf, int32_t *max_iter,
                          double *final_energy, int32_t *basis_set_n, int32_t *z_num);

void run_scf_c(int natoms, int nelectrons, int nalpha, int nbeta, const char *species, double *coords,
               double eps_scf, int max_iter, double *final_energy, int32_t *basis_set_n, int32_t *z_num) {
    run_scf_interface_py(&natoms, &nelectrons, &nalpha, &nbeta, species, coords,
                         &eps_scf, &max_iter, final_energy, basis_set_n, z_num);
}

