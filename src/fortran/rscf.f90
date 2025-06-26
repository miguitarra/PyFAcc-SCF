module rscf
    use constants_struct
    use molecule_struct
    use molecule
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
    use libint_integrals
    use iso_c_binding
    USE libint_f, ONLY: libint_t, libint2_static_init, libint2_static_cleanup, libint2_build, libint2_max_am_eri, &
                        compute_eri_f, libint2_init_eri, libint2_cleanup_eri
    implicit none

    private
    public :: run_rscf

contains

    subroutine run_rscf(mol, shells, basis_functions, eps_scf, max_iter, final_energy, basis_set_n)
        type(Molecule_type), intent(in) :: mol
        type(Shell_type), intent(in) :: shells(:)
        type(BasisFunction_type), intent(in) :: basis_functions(:)
        real(wp), intent(in)         :: eps_scf
        integer, intent(in)          :: max_iter, basis_set_n
        real(wp), intent(out)        :: final_energy
        type(ERI_Index_type), allocatable :: eri_list(:)

        integer :: mu, nu, lam, sig, nbf, iter, n, j, k, l, x, n_occ
        integer*8 :: n_unique, idx, i, nbf_sym
        real(wp) :: E_old, E_elec, delta_E, err_D, tei, sum_g, dummy_s   
        real(wp), allocatable :: S(:,:), T(:,:), V(:,:), H(:,:), S0(:,:)
        real(wp), allocatable :: D(:,:), F(:,:), G(:,:)
        real(wp), allocatable :: D_old(:,:)
        real(wp), allocatable :: eps(:)
        real(wp), allocatable :: C(:,:)
        real(wp), allocatable :: eigvals(:), eigvecs(:,:)
        real(wp), allocatable :: ERI(:,:,:,:), ERI_flat(:)
        integer, allocatable  :: index_map(:,:,:,:)
        real(wp) :: E_temp, ee_val1, ee_val2
        integer :: max_am

        max_am = shells(1)%l
        do i = 2, size(shells)
           if (shells(i)%l > max_am) max_am = shells(i)%l
        end do
        
        n = mol%nalpha + mol%nbeta
        n_occ = (n+1) / 2
                
        nbf = size(basis_functions)

        allocate(S(nbf,nbf), T(nbf,nbf), V(nbf,nbf), H(nbf,nbf))
        allocate(D(nbf,nbf), D_old(nbf,nbf))
        allocate(F(nbf,nbf), G(nbf,nbf))
        allocate(eps(nbf), C(nbf,nbf))
        allocate(eigvals(nbf))
        allocate(index_map(nbf, nbf, nbf, nbf))
        
        ! Build ERI index list
        nbf_sym = nbf*(nbf+1)/2
        n_unique = nbf_sym*(nbf_sym+1)/2

        allocate(eri_list(n_unique), ERI_flat(n_unique + 1))
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

        !acc data copyin(basis_functions, mol, mol%coords, mol%atomic_numbers) &
        !$acc data create(S, T, V, H, D, F, G, ERI_flat, index_map)

        ! One-electron integrals
        
        call S_overlap(nbf, basis_set_n , basis_functions, S)
        call T_kinetic(nbf, basis_set_n, basis_functions, T)
        call V_nuclear(nbf, basis_set_n, basis_functions, mol%coords, mol%atomic_numbers, mol%natoms, V)
        !print *, "One electron integrals calculated"
        

        !$acc wait

        H = T + V  ! Compute core matrix

        ! Compute ERIs

        call libint2_static_init()
        call compute_ERI_libint2(mol, shells,  basis_functions, n_unique, ERI_flat, eri_list, basis_set_n, index_map)
        !print *, ERI_flat
        !$acc update device(index_map, ERI_flat)

        ! Initial Guess
        call orthonormal_diagonalize(H, S, eigvals, C)
        D = 2*matmul(C(:,1:n_occ), transpose(C(:,1:n_occ)))

        E_old = 0.0_wp
        
        do iter = 1, max_iter
            G = 0.0_wp
            !$acc update device(D)
            !$acc parallel loop collapse(2) gang vector private(sum_g, ee_val1, ee_val2) present(D, G, ERI_flat, index_map)
            do mu = 1, nbf
                do nu = 1, nbf
                    sum_g = 0.0_wp
                    do lam = 1, nbf
                        do sig = 1, nbf
                            ee_val1 = get_ERI(mu, nu, lam, sig, ERI_flat, index_map)
                            ee_val2 = get_ERI(mu, sig, lam, nu, ERI_flat, index_map)
                            sum_g = sum_g + D(lam,sig) * (ee_val1 - 0.5_wp * ee_val2)
                        end do
                    end do
                    G(mu,nu) = sum_g
                end do
            end do
            !$acc update self(G)
        
            F = H + G
            D_old = D
            call orthonormal_diagonalize(F, S, eps, C)
            D = 2*matmul(C(:,1:n_occ), transpose(C(:,1:n_occ)))
        
            E_temp = 0.0_wp
            !$acc update device(D, H, F)
            !$acc parallel loop collapse(2) reduction(+:E_temp) present(D, H, F)
            do mu = 1, nbf
                do nu = 1, nbf
                    E_temp = E_temp + D(mu,nu) * (H(mu,nu) + F(mu,nu))
                end do
            end do
        
            E_elec = 0.5_wp * E_temp
            delta_E = abs(E_elec - E_old)
            err_D = sqrt(sum((D - D_old)**2))
        
            if (delta_E < eps_scf .and. err_D < eps_scf) exit
            E_old = E_elec
        end do


        
        !$acc exit data delete(basis_functions, mol, mol%coords, mol%atomic_numbers, S, T, V, H, D, F, G, ERI_flat, index_map)
        !$acc end data

        final_energy = E_elec + mol%Enuc

        deallocate(S, T, V, H, D, F, G, D_old, eps, C, ERI_flat, eri_list, eigvals, index_map)

    end subroutine run_rscf

    function get_ERI(mu, nu, lam, sig, ERI_flat, index_map) result(tei)
        !$acc routine seq
        integer, intent(in) :: mu, nu, lam, sig
        real(wp), dimension(:), intent(in) :: ERI_flat
        integer, intent(in) :: index_map(:,:,:,:)
        real(wp) :: tei
        integer :: idx
    
        idx = index_map(mu, nu, lam, sig)
        if (idx > 0) then
            tei = ERI_flat(idx)
        else
            tei = 0.0_wp
        end if
    end function get_ERI

end module rscf