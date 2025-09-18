#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// Fortran subroutine interface
void run_scf_interface_py(int32_t *natoms, int32_t *nelectrons, int32_t *charge, int32_t *nalpha, int32_t *nbeta,
                          const char *species, double *coords, double *eps_scf, int32_t *max_iter,
                          double *final_energy, const char *basis_set_c, int32_t *len_basis_set, int32_t *z_num);

void run_scf_c(int natoms, int nelectrons, int charge, int nalpha, int nbeta, const char *species, double *coords,
               double eps_scf, int max_iter, double *final_energy, const char *basis_set_c, int32_t *z_num) {
    int len_basis_set = strlen(basis_set_c);
    run_scf_interface_py(&natoms, &nelectrons, &charge, &nalpha, &nbeta, species, coords,
                         &eps_scf, &max_iter, final_energy, basis_set_c, &len_basis_set, z_num);
}

