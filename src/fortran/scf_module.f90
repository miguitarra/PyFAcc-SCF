module scf_module
    use integrals
    use constants_struct
    use molecule_module
    use molecule_struct
    use slater
    use utils_module
    use openacc
    implicit none

    type :: ERI_Index
        integer :: mu, nu, lam, sig
    end type ERI_Index

    private
    public :: run_scf

contains

    subroutine run_scf(mol, basis_functions, eps_scf, max_iter, final_energy)
        type(Molecule), intent(in) :: mol
        type(BasisFunction), intent(in) :: basis_functions(:)
        real(wp), intent(in)         :: eps_scf
        integer, intent(in)          :: max_iter
        real(wp), intent(out)        :: final_energy
        type(ERI_Index), allocatable :: eri_list(:)

        integer :: mu, nu, lam, sig, nbf, iter, n_alpha, n_beta, i, j, k, l, n_unique, idx
        real(wp) :: E_old, E_elec, delta_E, err_D, tei, sum_ga, sum_gb
        real(wp), allocatable :: S(:,:), T(:,:), V(:,:), H(:,:)
        real(wp), allocatable :: D_alpha(:,:), D_beta(:,:), F_alpha(:,:), F_beta(:,:), G_alpha(:,:), G_beta(:,:)
        real(wp), allocatable :: D_alpha_old(:,:), D_beta_old(:,:)
        real(wp), allocatable :: eps_alpha(:), eps_beta(:)
        real(wp), allocatable :: C_alpha(:,:), C_beta(:,:)
        real(wp), allocatable :: eigvals(:), eigvecs(:,:), centers(:,:)
        real(wp), allocatable :: ERI(:,:,:,:), ERI_flat(:)
        real(wp) :: E_temp

        n_alpha = mol%nalpha
        n_beta  = mol%nbeta
        nbf = size(basis_functions)

        allocate(S(nbf,nbf), T(nbf,nbf), V(nbf,nbf), H(nbf,nbf))
        allocate(D_alpha(nbf,nbf), D_beta(nbf,nbf), D_alpha_old(nbf,nbf), D_beta_old(nbf,nbf))
        allocate(F_alpha(nbf,nbf), F_beta(nbf,nbf), G_alpha(nbf,nbf), G_beta(nbf,nbf))
        allocate(eps_alpha(nbf), eps_beta(nbf), C_alpha(nbf,nbf), C_beta(nbf,nbf))
        allocate(eigvals(nbf), eigvecs(nbf,nbf), centers(3, nbf))
        allocate(ERI(nbf,nbf,nbf,nbf))

        ! Build ERI index list
        n_unique = 0
        do i = 1, nbf
            do j = 1, i
                do k = 1, nbf
                    do l = 1, k
                        if ((i*(i-1)/2 + j) >= (k*(k-1)/2 + l)) then
                            n_unique = n_unique + 1
                        end if
                    end do
                end do
            end do
        end do

        allocate(eri_list(n_unique), ERI_flat(n_unique))
        ERI_flat = 0.0_wp

        idx = 0
        do i = 1, nbf
            do j = 1, i
                do k = 1, nbf
                    do l = 1, k
                        if ((i*(i-1)/2 + j) >= (k*(k-1)/2 + l)) then
                            idx = idx + 1
                            eri_list(idx)%mu  = i
                            eri_list(idx)%nu  = j
                            eri_list(idx)%lam = k
                            eri_list(idx)%sig = l
                        end if
                    end do
                end do
            end do
        end do

        !$acc data copyin(basis_functions, mol, mol%coords, mol%atomic_numbers, eri_list) &
        !$acc     create(S, T, V, H, D_alpha, D_beta, F_alpha, F_beta, G_alpha, G_beta, ERI, ERI_flat)
        do i = 1, size(basis_functions)
            !$acc enter data copyin(basis_functions(i)%exponents, basis_functions(i)%coefficients, basis_functions(i)%center)
        end do

            ! One-electron integrals
            !$acc parallel loop collapse(2) present(mol, basis_functions, S, T, V)
            do mu = 1, nbf
                do nu = 1, nbf
                    call oneint(mol%coords, real(mol%atomic_numbers, wp), &
                                 basis_functions(mu)%center, basis_functions(nu)%center, &
                                 basis_functions(mu)%exponents, basis_functions(nu)%exponents, &
                                 basis_functions(mu)%coefficients, basis_functions(nu)%coefficients, &
                                 S(mu,nu), T(mu,nu), V(mu,nu))
                end do
            end do
            !$acc update self(S, T, V)

            H = T + V

            ! Compute ERIs
            !$acc parallel loop present(eri_list, mol, basis_functions, ERI_flat)
            do i = 1, n_unique
                mu  = eri_list(i)%mu
                nu  = eri_list(i)%nu
                lam = eri_list(i)%lam
                sig = eri_list(i)%sig

                call twoint(mol%coords(:,basis_functions(mu)%atom), mol%coords(:,basis_functions(nu)%atom), &
                            mol%coords(:,basis_functions(lam)%atom), mol%coords(:,basis_functions(sig)%atom), &
                            basis_functions(mu)%exponents, basis_functions(nu)%exponents, &
                            basis_functions(lam)%exponents, basis_functions(sig)%exponents, &
                            basis_functions(mu)%coefficients, basis_functions(nu)%coefficients, &
                            basis_functions(lam)%coefficients, basis_functions(sig)%coefficients, &
                            tei)

                ERI_flat(i) = tei
            end do

            !$acc parallel loop present(ERI, ERI_flat, eri_list)
            do i = 1, n_unique
                mu  = eri_list(i)%mu
                nu  = eri_list(i)%nu
                lam = eri_list(i)%lam
                sig = eri_list(i)%sig
                tei = ERI_flat(i)

                ERI(mu,nu,lam,sig) = tei
                ERI(nu,mu,lam,sig) = tei
                ERI(mu,nu,sig,lam) = tei
                ERI(nu,mu,sig,lam) = tei
                ERI(lam,sig,mu,nu) = tei
                ERI(sig,lam,mu,nu) = tei
                ERI(lam,sig,nu,mu) = tei
                ERI(sig,lam,nu,mu) = tei
            end do
            !$acc update self(ERI)

            ! Initial Guess
            call orthonormal_diagonalize(H, S, eigvals, eigvecs)
            C_alpha = eigvecs
            C_beta  = eigvecs
            D_alpha = matmul(C_alpha(:,1:n_alpha), transpose(C_alpha(:,1:n_alpha)))
            D_beta  = matmul(C_beta(:,1:n_beta),  transpose(C_beta(:,1:n_beta)))

            E_old = 0.0_wp
            do iter = 1, max_iter
                G_alpha = 0.0_wp
                G_beta  = 0.0_wp

                !$acc update device(D_alpha, D_beta)
                !$acc parallel loop collapse(2) private(sum_ga, sum_gb) present(D_alpha, D_beta, ERI, G_alpha, G_beta)
                do mu = 1, nbf
                    do nu = 1, nbf
                        sum_ga = 0.0_wp
                        sum_gb = 0.0_wp
                        do lam = 1, nbf
                            do sig = 1, nbf
                                sum_ga = sum_ga + (D_alpha(lam,sig) + D_beta(lam,sig)) * ERI(mu,nu,lam,sig) - &
                                                D_alpha(lam,sig) * ERI(mu,sig,lam,nu)
                                sum_gb = sum_gb + (D_alpha(lam,sig) + D_beta(lam,sig)) * ERI(mu,nu,lam,sig) - &
                                                D_beta(lam,sig) * ERI(mu,sig,lam,nu)
                            end do
                        end do
                        G_alpha(mu,nu) = sum_ga
                        G_beta(mu,nu)  = sum_gb
                    end do
                end do
                !$acc update self(G_alpha, G_beta)

                F_alpha = H + G_alpha
                F_beta  = H + G_beta

                D_alpha_old = D_alpha
                D_beta_old = D_beta

                call orthonormal_diagonalize(F_alpha, S, eps_alpha, C_alpha)
                call orthonormal_diagonalize(F_beta,  S, eps_beta,  C_beta)

                D_alpha = matmul(C_alpha(:,1:n_alpha), transpose(C_alpha(:,1:n_alpha)))
                D_beta  = matmul(C_beta(:,1:n_beta),  transpose(C_beta(:,1:n_beta)))

                E_temp = 0.0_wp
                !$acc update device(D_alpha, D_beta, H, F_alpha, F_beta)
                !$acc parallel loop collapse(2) reduction(+:E_temp) present(D_alpha, D_beta, H, F_alpha, F_beta)
                do mu = 1, nbf
                    do nu = 1, nbf
                        E_temp = E_temp + 0.5_wp * (D_alpha(mu,nu) + D_beta(mu,nu)) * H(mu,nu)
                        E_temp = E_temp + 0.5_wp * (D_alpha(mu,nu) * F_alpha(mu,nu) + D_beta(mu,nu) * F_beta(mu,nu))
                    end do
                end do
                E_elec = E_temp

                delta_E = abs(E_elec - E_old)
                err_D = sqrt(sum((D_alpha - D_alpha_old)**2 + (D_beta - D_beta_old)**2))

                if (delta_E < eps_scf .and. err_D < eps_scf) exit
                E_old = E_elec
            end do
        do i = 1, size(basis_functions)
            !$acc exit data delete(basis_functions(i)%exponents, basis_functions(i)%coefficients, basis_functions(i)%center)
        end do
        !$acc exit data delete(basis_functions, mol, mol%coords, mol%atomic_numbers, eri_list,S, T, V, H, D_alpha, D_beta, F_alpha, F_beta, G_alpha, G_beta, ERI, ERI_flat)
        !$acc end data

        final_energy = E_elec + mol%Enuc

        deallocate(S, T, V, H, D_alpha, D_beta, F_alpha, F_beta, G_alpha, G_beta, &
                   D_alpha_old, D_beta_old, eps_alpha, eps_beta, C_alpha, C_beta, ERI, ERI_flat, eri_list)

    end subroutine run_scf

    subroutine orthonormal_diagonalize(F, S, eigvals, eigvecs)
        real(wp), intent(in)  :: F(:,:), S(:,:)
        real(wp), intent(out) :: eigvals(:), eigvecs(:,:)
        real(wp), allocatable :: S_copy(:,:), S_eigvals(:), S_eigvecs(:,:), X(:,:), F_prime(:,:)
        integer :: n, i, j

        n = size(F,1)
        allocate(S_copy(n,n), S_eigvals(n), S_eigvecs(n,n), X(n,n), F_prime(n,n))

        S_copy = S
        call diagonalize_symmetric(S_copy, S_eigvals, S_eigvecs)

        do j = 1, n
            if (S_eigvals(j) < 1.0e-8_wp) then
                print *, "Warning: Small overlap eigenvalue detected: ", S_eigvals(j)
            end if
        end do

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

end module scf_module
