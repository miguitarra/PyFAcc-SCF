module scf_module
    use integrals
    use molecule_module
    use slater
    use utils_module
    implicit none
    private
    public :: run_scf

contains 

    subroutine run_scf(mol, basis_functions, eps_scf, max_iter, final_energy)
        type(Molecule), intent(in) :: mol
        type(BasisFunction), intent(in) :: basis_functions(:)
        real(wp), intent(in)         :: eps_scf
        integer, intent(in)          :: max_iter
        real(wp), intent(out)        :: final_energy

        integer :: mu, nu, lam, sig, nbf, iter, n_occ
        real(wp) :: E_old, E_elec, delta_E, err_D, tei
        real(wp), allocatable :: S(:,:), T(:,:), V(:,:), H(:,:), D(:,:), F(:,:), G(:,:), D_old(:,:), eps(:), C(:,:)
        real(wp), allocatable :: eigvals(:), eigvecs(:,:)
        real(wp), allocatable :: centers(:,:)
        real(wp), allocatable :: ERI(:,:,:,:)
        logical :: is_closed_shell
        is_closed_shell = (mod(mol%nelectrons, 2) == 0)

        ! Setup dimensions
        nbf = size(basis_functions)

        allocate(S(nbf,nbf), T(nbf,nbf), V(nbf,nbf), H(nbf,nbf), D(nbf,nbf), F(nbf,nbf), G(nbf,nbf))
        allocate(D_old(nbf,nbf), eps(nbf), C(nbf,nbf))
        allocate(eigvals(nbf), eigvecs(nbf,nbf), centers(3, nbf))
        allocate(ERI(nbf,nbf,nbf,nbf))

        n_occ = (mol%nelectrons + 1) / 2

        ! One electron integrals
        do mu = 1, nbf
            do nu = 1, nbf
                !print *, "Computing integrals for i = ", mu, ", j = ", nu
                call oneint( &
                    mol%coords, &
                    real(mol%atomic_numbers, wp), &
                    basis_functions(mu)%center, &
                    basis_functions(nu)%center, &
                    basis_functions(mu)%exponents, &
                    basis_functions(nu)%exponents, &
                    basis_functions(mu)%coefficients, &
                    basis_functions(nu)%coefficients, &
            S(mu,nu), T(mu,nu), V(mu,nu))
            end do
        end do

        !print *, "Sum of S matrix elements:", sum(S)
        !print *, "Sum of T matrix elements:", sum(T)
        !print *, "Sum of V matrix elements:", sum(V)

        H = T + V
        !print *, "One-electron integrals built."

        ! Two-electron integrals
        ERI = 0.0_wp
        do mu = 1, nbf
            do nu = 1, nbf
                do lam = 1, nbf
                    do sig = 1, nbf
                        call twoint( &
                            mol%coords(:,basis_functions(mu)%atom), &
                            mol%coords(:,basis_functions(nu)%atom), &
                            mol%coords(:,basis_functions(lam)%atom), &
                            mol%coords(:,basis_functions(sig)%atom), &
                            basis_functions(mu)%exponents, &
                            basis_functions(nu)%exponents, &
                            basis_functions(lam)%exponents, &
                            basis_functions(sig)%exponents, &
                            basis_functions(mu)%coefficients, &
                            basis_functions(nu)%coefficients, &
                            basis_functions(lam)%coefficients, &
                            basis_functions(sig)%coefficients, &
                            tei)
                        ERI(mu,nu,lam,sig) = tei
                    end do
                end do
            end do
        end do
        ! print *, "ERI(1,1,1,1):", ERI(1,1,1,1)
        !print *, "Two-electron integrals built."

        ! Initial guess
        call orthonormal_diagonalize(H, S, eps, C)
        !print *, "Initial eigenvalues: ", eps
        if (is_closed_shell) then ! will have to change this for open shell
            D = 2*matmul(C(:,1:n_occ), transpose(C(:,1:n_occ)))
        else 
            D = matmul(C(:,1:n_occ), transpose(C(:,1:n_occ)))
        end if
        !print *, "Initial guess density matrix built."

    
        E_old = 0.0_wp
        do iter = 1, max_iter
            ! Build G matrix
            if (mol%nelectrons > 1) then
                ! Build G matrix (only for multi-electron systems)
                G = 0.0_wp
                do mu = 1, nbf
                    do nu = 1, nbf
                        do lam = 1, nbf
                            do sig = 1, nbf
                                G(mu,nu) = G(mu,nu) + D(lam,sig) * (ERI(mu,nu,lam,sig) - 0.5_wp * ERI(mu,sig,lam,nu))
                            end do
                        end do
                    end do
                end do
            else
                G = 0.0_wp  ! no self-interaction
            end if
            F = H + G

            ! Orthonormalize and Diagonalize
            call orthonormal_diagonalize(F, S, eigvals, eigvecs)

            ! Build new density matrix
            D_old = D
            C = eigvecs
            if (is_closed_shell) then ! will have to change this for open shell
                D = 2*matmul(C(:,1:n_occ), transpose(C(:,1:n_occ)))
            else 
                D = matmul(C(:,1:n_occ), transpose(C(:,1:n_occ)))
            end if
    
            ! Compute energy
            E_elec = 0.0_wp
            do mu = 1, nbf
                do nu = 1, nbf
                    E_elec = E_elec + D(mu,nu)*(H(mu,nu) + F(mu,nu))
                end do
            end do

            !Convergence
            delta_E = abs(E_elec - E_old)
            err_D = sqrt(sum((D - D_old)**2))
            !write(*,'(A,I6,A,F16.10,A,E10.2E2)') "Iter ", iter, "  E = ", E_elec/2 + mol%Enuc , "  dE = ", delta_E  ! /2 term
            if (delta_E < eps_scf) exit
            E_old = E_elec
        end do
        final_energy = E_elec/2 + mol%Enuc ! /2 term

        !write(*,'(A,F16.10)') "Final SCF Energy: ", final_energy 

        deallocate(S, T, V, H, D, F, G, eigvals, eigvecs, centers)

    contains 
        subroutine orthonormal_diagonalize(F, S, eigvals, eigvecs)
            real(wp), intent(in)  :: F(:,:), S(:,:)
            real(wp), intent(out) :: eigvals(:), eigvecs(:,:)
            real(c_double), allocatable :: S_copy(:,:), S_eigvals(:), S_eigvecs(:,:), X(:,:), F_prime(:,:)
            integer :: n, i, j

            n = size(F,1)
            allocate(S_copy(n,n), S_eigvals(n), S_eigvecs(n,n), X(n,n), F_prime(n,n))

            S_copy = S
            call diagonalize_symmetric(S_copy, S_eigvals, S_eigvecs)

            do i = 1, n
                do j = 1, n
                    X(i,j) = S_eigvecs(i,j) / sqrt(S_eigvals(j))
                end do
            end do

            F_prime = matmul(transpose(X), matmul(F, X))
            call diagonalize_symmetric(F_prime, eigvals, eigvecs)
            eigvecs = matmul(X, eigvecs)

            deallocate(S_copy, S_eigvals, S_eigvecs, X, F_prime)
        end subroutine orthonormal_diagonalize

end subroutine run_scf

end module scf_module