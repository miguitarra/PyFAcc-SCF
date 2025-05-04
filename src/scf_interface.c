#include <stdio.h>

extern void scf_h2(int *natoms, int *atomic_numbers, double *coords, double *eps_scf, int *max_iter, double *final_energy, int *nbasis_atom, double *basis_exp);


void run_scf_interface_py(int natoms, int *atomic_numbers, double *coords, double eps_scf, int max_iter, double *final_energy, int *nbasis_atom, double *basis_exp) {
    // Call Fortran subroutine run_scf()
    scf_h2(&natoms, atomic_numbers, coords, &eps_scf, &max_iter, final_energy, nbasis_atom, basis_exp);
}
