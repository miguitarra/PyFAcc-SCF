#pragma once
#include <cstddef>

extern "C" {

/**
 * Compute one-electron integrals using Libint
 * @param n_shells Number of shells
 * @param n_basis Number of basis functions (total)
 * @param centers 3*n_shells array of shell centers (row major)
 * @param exponents Variable length, flattened exponents (see shell_offsets)
 * @param coefficients Variable length, flattened coefficients (see shell_offsets)
 * @param l_shells Angular momentum of each shell (size n_shells)
 * @param shell_offsets Start index into exponents/coefficients for each shell
 * @param n_primitives Number of primitives per shell (size n_shells)
 * @param result Output: the computed integrals (size n_basis*n_basis)
 */
void libint_calc_1body_ints(
    int n_shells,
    int n_basis,
    const double* centers,
    const double* exponents,
    const double* coefficients,
    const int* l_shells,
    const int* shell_offsets,
    const int* n_primitives,
    double* result,
    int op_type  // 0=overlap, 1=kinetic, 2=nuclear, etc.
);

}