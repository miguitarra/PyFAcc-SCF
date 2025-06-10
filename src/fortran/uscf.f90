module uscf
    use constants_struct
    use molecule_struct
    use molecule_struct
    use basis_function_struct
    use utils
    use openacc
    use la
    use io
    use nuclear
    use overlap
    use kinetic
    use eri_struct
    use electronic
    use integrals
    implicit none


    private
    public :: run_uscf


contains

    subroutine run_uscf(mol, shells, basis_functions, eps_scf, max_iter, final_energy, basis_set_n)
        type(Molecule_type), intent(in) :: mol
        type(Shell_type), intent(in) :: shells(:)
        type(BasisFunction_type), intent(in) :: basis_functions(:)
        real(wp), intent(in)         :: eps_scf
        integer, intent(in)          :: max_iter, basis_set_n
        real(wp), intent(out)        :: final_energy
        type(ERI_Index_type), allocatable :: eri_list(:)
        

        integer :: mu, nu, lam, sig, nbf, iter, n_alpha, n_beta, j, k, l, x
        integer*8 :: n_unique, idx, i, nbf_sym
        real(wp) :: E_old, E_elec, delta_E, err_D, tei, sum_ga, sum_gb, dummy_s   
        real(wp), allocatable :: S(:,:), T(:,:), V(:,:), H(:,:), S0(:,:)
        real(wp), allocatable :: D_alpha(:,:), D_beta(:,:), F_alpha(:,:), F_beta(:,:), G_alpha(:,:), G_beta(:,:)
        real(wp), allocatable :: D_alpha_old(:,:), D_beta_old(:,:)
        real(wp), allocatable :: eps_alpha(:), eps_beta(:)
        real(wp), allocatable :: C_alpha(:,:), C_beta(:,:)
        real(wp), allocatable :: eigvals(:), eigvecs(:,:)
        real(wp), allocatable :: ERI_flat(:)
        integer, allocatable :: eri_map(:,:)
        real(wp) :: E_temp, ee_val1, ee_val2


        ! --- Check where it is being run ---
        integer :: dev_type
        dev_type = acc_get_device_type()

        select case (dev_type)
        case (acc_device_nvidia)
            !print *, "OpenACC: Running on NVIDIA GPU."
        case (acc_device_host)
            !Only let know if it is running on host
            print *, "OpenACC: Running on host (CPU)."
        case default
            !print *, "OpenACC: Unknown or unsupported device type."
        end select
        ! ---------------------------

        n_alpha = mol%nalpha
        n_beta  = mol%nbeta
        nbf = size(basis_functions)

        allocate(S(nbf,nbf), T(nbf,nbf), V(nbf,nbf), H(nbf,nbf))
        allocate(D_alpha(nbf,nbf), D_beta(nbf,nbf), D_alpha_old(nbf,nbf), D_beta_old(nbf,nbf))
        allocate(F_alpha(nbf,nbf), F_beta(nbf,nbf), G_alpha(nbf,nbf), G_beta(nbf,nbf))
        allocate(eps_alpha(nbf), eps_beta(nbf), C_alpha(nbf,nbf), C_beta(nbf,nbf))
        allocate(eigvals(nbf), eigvecs(nbf,nbf))

        !do x = 1, nbf
        !    print *, "Basis function:", x
        !    print *, "Ang mom"
        !    print *, basis_functions(x)%ang_mom
        !    print *, "Exponents"
        !    call print_vector(basis_functions(x)%exponents, 6)
        !    print *, "Coefficients"
        !    call print_vector(basis_functions(x)%coefficients, 6)
        !    print *, "Centers"
        !    call print_vector( basis_functions(x)%center, 3)
        !end do

        
        nbf_sym = nbf*(nbf+1)/2
        n_unique = nbf_sym*(nbf_sym+1)/2

        print *, "Total unique ERIs: ", n_unique

        ! Second pass: populate eri_map (4 x n_unique)
        allocate(eri_map(4, n_unique))
        idx = 0
        !acc parallel loop collapse(4) present(eri_map) private(idx)
        do i = 1, nbf
            do j = 1, i
                do k = 1, nbf
                    do l = 1, k
                        if ((i*(i-1)/2 + j) >= (k*(k-1)/2 + l)) then
                            !acc atomic capture
                            idx = idx + 1
                            eri_map(1, idx) = i
                            eri_map(2, idx) = j
                            eri_map(3, idx) = k
                            eri_map(4, idx) = l
                        end if
                    end do
                end do
            end do
        end do

        ! Now populate eri_list from eri_map
        allocate(eri_list(n_unique), ERI_flat(n_unique+1))
        ERI_flat = 0.0_wp

        !acc parallel loop present(eri_list, eri_map)
        do i = 1, n_unique
            eri_list(i)%mu  = eri_map(1, i)
            eri_list(i)%nu  = eri_map(2, i)
            eri_list(i)%lam = eri_map(3, i)
            eri_list(i)%sig = eri_map(4, i)
        end do
        deallocate(eri_map)
        

        print *, "ERI list populated"

        !acc data copyin(basis_functions, mol, mol%atomic_numbers, eri_list) &
        !acc     create(S, T, V, H, D_alpha, D_beta, F_alpha, F_beta, G_alpha, G_beta, ERI_flat)
        do i = 1, size(basis_functions)
            !acc enter data copyin(basis_functions(i)%exponents, basis_functions(i)%coefficients, basis_functions(i)%center, basis_functions(i)%ang_mom)
        end do

        ! One-electron integrals
        call S_overlap(nbf, basis_set_n , basis_functions, S)
        call T_kinetic(nbf, basis_set_n, basis_functions, T)
        call V_nuclear(nbf, basis_set_n, basis_functions, mol%coords, mol%atomic_numbers, mol%natoms, V)

        print *, "Full S matrix:" ! Overlap
        call print_matrix(S, nbf)
        print *, "Full T matrix:" ! Kinetic 
        call print_matrix(T, nbf)
        print *, "Full V matrix:" ! Potential
        call print_matrix(V, nbf)

        H = T + V  ! Compute core matrix
        print *, "One electron integrals computed"

        ! Compute ERIs
        call compute_ERI_libint(mol, shells,  basis_functions, n_unique, ERI_flat, eri_list, basis_set_n)
        print *, "Two-electron integrals computed"

        ! Initial Guess
        call orthonormal_diagonalize(H, S, eigvals, eigvecs)
        C_alpha = eigvecs
        C_beta  = eigvecs
        D_alpha = matmul(C_alpha(:,1:n_alpha), transpose(C_alpha(:,1:n_alpha)))
        D_beta  = matmul(C_beta(:,1:n_beta),  transpose(C_beta(:,1:n_beta)))

        E_old = 0.0_wp
        print *, "Entering scf loop"
        do iter = 1, max_iter
            G_alpha = 0.0_wp
            G_beta  = 0.0_wp

            !acc update device(D_alpha, D_beta)
            !acc parallel loop collapse(2) private(sum_ga, sum_gb) present(D_alpha, D_beta, G_alpha, G_beta)
            do mu = 1, nbf
                do nu = 1, nbf
                    sum_ga = 0.0_wp
                    sum_gb = 0.0_wp
                    do lam = 1, nbf
                        do sig = 1, nbf
                            ee_val1 = get_ERI(mu, nu, lam, sig, ERI_flat)
                            ee_val2 = get_ERI(mu, sig, lam, nu, ERI_flat)
                            sum_ga = sum_ga + (D_alpha(lam,sig) + D_beta(lam,sig)) * ee_val1 - &
                                            D_alpha(lam,sig) * ee_val2
                            sum_gb = sum_gb + (D_alpha(lam,sig) + D_beta(lam,sig)) * ee_val1 - &
                                            D_beta(lam,sig) * ee_val2
                        end do
                    end do
                    G_alpha(mu,nu) = sum_ga
                    G_beta(mu,nu)  = sum_gb
                end do
            end do
            !acc update self(G_alpha, G_beta)

            F_alpha = H + G_alpha
            F_beta  = H + G_beta

            D_alpha_old = D_alpha
            D_beta_old = D_beta

            call orthonormal_diagonalize(F_alpha, S, eps_alpha, C_alpha)
            call orthonormal_diagonalize(F_beta,  S, eps_beta,  C_beta)

            D_alpha = matmul(C_alpha(:,1:n_alpha), transpose(C_alpha(:,1:n_alpha)))
            D_beta  = matmul(C_beta(:,1:n_beta),  transpose(C_beta(:,1:n_beta)))

            E_temp = 0.0_wp
            !acc update device(D_alpha, D_beta, H, F_alpha, F_beta)
            !acc parallel loop collapse(2) reduction(+:E_temp) present(D_alpha, D_beta, H, F_alpha, F_beta) 
            do mu = 1, nbf
                do nu = 1, nbf
                    E_temp = E_temp + 0.5_wp * (D_alpha(mu,nu) + D_beta(mu,nu)) * H(mu,nu)
                    E_temp = E_temp + 0.5_wp * (D_alpha(mu,nu) * F_alpha(mu,nu) + D_beta(mu,nu) * F_beta(mu,nu))
                end do
            end do
            E_elec = E_temp

            delta_E = abs(E_elec - E_old)
            err_D = sqrt(sum((D_alpha - D_alpha_old)**2 + (D_beta - D_beta_old)**2))
            write(*,'(A,I6,A,F16.10,A,E10.2E2)') "Iter ", iter, "  E = ", E_elec + mol%Enuc, "  dE = ", delta_E


            
            if (delta_E < eps_scf .and. err_D < eps_scf) exit
            
            
            E_old = E_elec
        end do
        
        do i = 1, size(basis_functions)
            !acc exit data delete(basis_functions(i)%exponents, basis_functions(i)%coefficients, basis_functions(i)%center)
        end do
        !acc exit data delete(basis_functions, mol, mol%atomic_numbers, eri_list,S, T, V, H, D_alpha, D_beta, F_alpha, F_beta, G_alpha, G_beta, ERI_flat)
        !acc end data

        final_energy = E_elec + mol%Enuc


        deallocate(S, T, V, H, D_alpha, D_beta, F_alpha, F_beta, G_alpha, G_beta, &
                   D_alpha_old, D_beta_old, eps_alpha, eps_beta, C_alpha, C_beta, ERI_flat, eri_list)

    end subroutine run_uscf

    function get_ERI(mu, nu, lam, sig, ERI_flat) result(tei)
        integer, intent(in) :: mu, nu, lam, sig
        real(wp), dimension(:), intent(in) :: ERI_flat
        real(wp) :: tei
        integer :: idx
    
        idx = index_map(mu, nu, lam, sig)
        if (idx > 0) then
            tei = ERI_flat(idx)
        else
            tei = 0.0_wp
        end if
    end function get_ERI





end module uscf