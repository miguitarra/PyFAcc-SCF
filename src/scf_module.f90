module scf_module
    use molecule_module
    use integral_module
    use utils_module
    use iso_c_binding
    implicit none

contains

    subroutine run_scf(mol, basis_functions, eps_scf, max_iter, final_energy)
        type(Molecule), intent(in) :: mol 
        type(BasisFunction), intent(in) :: basis_functions(:)
        real(c_double), intent(in) :: eps_scf
        integer, intent(in) :: max_iter
        real(c_double), intent(out) :: final_energy

        integer :: nbf, iter
        real(c_double) :: energy_old, energy_new, delta_energy
        real(c_double), allocatable :: S(:,:), Hcore(:,:), ERIs(:,:,:,:), F(:,:), C(:,:), D(:,:), eigvals(:)
        integer :: i, j, k, l

        nbf = size(basis_functions)

        call build_integrals(mol, basis_functions, S, Hcore, ERIs) ! Rewrite this to use the new module

        allocate(F(nbf,nbf), C(nbf,nbf), D(nbf,nbf), eigvals(nbf))
        D = 0.0d0
        F = Hcore

        iter = 0
        energy_old = 0.0d0

        do while (iter < max_iter)

            ! Diagonalize F
            call diagonalize_symmetric(F, eigvals, C)

            ! Build density matrix (only 1 occupied orbital)
            D = 2.0d0 * matmul(C(:,1:1), transpose(C(:,1:1)))

            ! Build new Fock matrix
            F = Hcore
            do i = 1, nbf
                do j = 1, nbf
                    do k = 1, nbf
                        do l = 1, nbf
                            F(i,j) = F(i,j) + D(k,l) * (2.0d0 * ERIs(i,k,j,l) - ERIs(i,k,l,j))
                        end do
                    end do
                end do
            end do

            ! Compute electronic energy
            energy_new = 0.0d0
            do i = 1, nbf
                do j = 1, nbf
                    energy_new = energy_new + 0.5d0 * D(i,j) * (Hcore(i,j) + F(i,j))
                end do
            end do

            delta_energy = abs(energy_new - energy_old)
            if (delta_energy < eps_scf) exit

            energy_old = energy_new
            iter = iter + 1

        end do

        final_energy = energy_new

    end subroutine run_scf

end module scf_module
