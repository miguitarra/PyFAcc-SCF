#include "libint_wrapper.h"
#include <libint2.hpp>
#include <vector>
#include <cassert>

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
    int op_type
) {
    using libint2::Shell;
    using libint2::Operator;
    
    std::vector<Shell> shells;
    for (int i = 0; i < n_shells; ++i) {
        std::vector<double> exps(exponents + shell_offsets[i], exponents + shell_offsets[i] + n_primitives[i]);
        std::vector<double> coefs(coefficients + shell_offsets[i], coefficients + shell_offsets[i] + n_primitives[i]);
        std::vector<double> center(centers + 3*i, centers + 3*i + 3);
        shells.emplace_back(l_shells[i], exps, {coefs}, center);
    }

    libint2::Engine engine(
        op_type == 0 ? Operator::overlap :
        op_type == 1 ? Operator::kinetic :
        op_type == 2 ? Operator::nuclear : Operator::overlap,  // Expand as needed
        shells.max_nprim(), shells.max_l()
    );

    // For nuclear, also set positions/charges here if needed

    std::fill(result, result + n_basis * n_basis, 0.0);
    int bf1 = 0;
    for (int s1 = 0; s1 < n_shells; ++s1) {
        int n1 = shells[s1].size();
        int bf2 = 0;
        for (int s2 = 0; s2 < n_shells; ++s2) {
            int n2 = shells[s2].size();
            const auto& buf = engine.compute(shells[s1], shells[s2]);
            if (buf != nullptr) {
                for (int f1 = 0; f1 < n1; ++f1)
                for (int f2 = 0; f2 < n2; ++f2)
                    result[(bf1 + f1) * n_basis + (bf2 + f2)] = buf[f1 * n2 + f2];
            }
            bf2 += n2;
        }
        bf1 += n1;
    }
}
