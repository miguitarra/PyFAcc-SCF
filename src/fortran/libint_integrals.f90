module libint_integrals
    use ISO_C_BINDING, ONLY: C_DOUBLE, C_F_POINTER, C_NULL_PTR
    USE libint_f, ONLY: libint_t, libint2_static_init, libint2_static_cleanup, libint2_build, libint2_max_am_eri, &
                        compute_eri_f, libint2_init_eri, libint2_cleanup_eri

    use molecule_struct
    use basis_function_struct
    use electronic
    use utils
    use boys
    use constants_struct

    implicit none

    integer, private :: index_map_nbf = 0
    real(C_DOUBLE), parameter :: ERI_MIN = -1.0d3
    real(C_DOUBLE), parameter :: ERI_MAX =  1.0d3

    interface
    
        !----------------------------------------------------------
        ! One electron Integrals Wrapper
        !----------------------------------------------------------
        subroutine libint_overlap(amA, contrdepthA, A, alphaA, cA, &
                                  amB, contrdepthB, B, alphaB, cB, &
                                  overlapAB) bind(C, name="libint_overlap_c")
          import :: C_INT, C_DOUBLE
          integer(C_INT), intent(in)  :: amA, contrdepthA
          real(C_DOUBLE), intent(in) :: A(*)
          real(C_DOUBLE), intent(in) :: alphaA(*)
          real(C_DOUBLE), intent(in) :: cA(*)
          integer(C_INT), intent(in)  :: amB, contrdepthB
          real(C_DOUBLE), intent(in) :: B(*)
          real(C_DOUBLE), intent(in) :: alphaB(*)
          real(C_DOUBLE), intent(in) :: cB(*)
          real(C_DOUBLE), intent(out) :: overlapAB(*)
    
        end subroutine libint_overlap
    
        subroutine libint_kinetic(amA, contrdepthA, A, alphaA, cA, &
                                  amB, contrdepthB, B, alphaB, cB, &
                                  kineticAB) bind(C, name="libint_kinetic_c")
          import :: C_INT, C_DOUBLE
          integer(C_INT), intent(in)  :: amA, contrdepthA
          real(C_DOUBLE), intent(in)  :: A(*)
          real(C_DOUBLE), intent(in)  :: alphaA(*)
          real(C_DOUBLE), intent(in)  :: cA(*)
          integer(C_INT), intent(in)  :: amB, contrdepthB
          real(C_DOUBLE), intent(in)  :: B(*)
          real(C_DOUBLE), intent(in)  :: alphaB(*)
          real(C_DOUBLE), intent(in)  :: cB(*)
          real(C_DOUBLE), intent(out) :: kineticAB(*)
    
        end subroutine libint_kinetic
    
        subroutine libint_potential(amA, contrdepthA, A, alphaA, cA, &
                                  amB, contrdepthB, B, alphaB, cB, &
                                  C, potAB) bind(C, name="libint_potential_c")
          import :: C_INT, C_DOUBLE
          integer(C_INT), intent(in)    :: amA, contrdepthA
          real(C_DOUBLE), intent(in)    :: A(*)
          real(C_DOUBLE), intent(in)    :: alphaA(*)
          real(C_DOUBLE), intent(in)    :: cA(*)
          integer(C_INT), intent(in)    :: amB, contrdepthB
          real(C_DOUBLE), intent(in)    :: B(*)
          real(C_DOUBLE), intent(in)    :: alphaB(*)
          real(C_DOUBLE), intent(in)    :: cB(*)
          real(C_DOUBLE), intent(in)    :: C(*)
          real(C_DOUBLE), intent(inout) :: potAB(*)
    
        end subroutine libint_potential



    end interface 

contains

    integer function nint_cartesian(am)
        implicit none
        integer, intent(in) :: am
        nint_cartesian = (am+1)*(am+2)/2
    end function nint_cartesian

    subroutine compute_S_libint(mol, shells, basis_function, basis_set_n, S)
        type(Molecule_type), intent(in) :: mol
        type(Shell_type), intent(in) :: shells(:)
        type(BasisFunction_type), intent(in) :: basis_function(:)
        integer, intent(in) :: basis_set_n
        real(wp), intent(out) :: S(:,:)

        integer :: s1, s2, i, j
        integer :: am1, am2, n1, n2
        integer :: ncart1, ncart2
        integer :: bf1_start, bf2_start
        integer :: idx1, idx2

        real(c_double) :: A(3), B(3)
        real(c_double) :: exp1(basis_set_n), exp2(basis_set_n)
        real(c_double) :: coeff1(basis_set_n), coeff2(basis_set_n)
        real(c_double), allocatable :: overlapAB(:)

        do s1 = 1, size(shells)
            am1 = shells(s1)%l
            A = shells(s1)%center
            exp1 = shells(s1)%exponents
            coeff1 = shells(s1)%coefficients
            ncart1 = nint_cartesian(am1)
    
            do s2 = s1, size(shells)
                am2 = shells(s2)%l
                B = shells(s2)%center
                exp2 = shells(s2)%exponents
                coeff2 = shells(s2)%coefficients
                ncart2 = nint_cartesian(am2)
    
                allocate(overlapAB(ncart1 * ncart2))

                print *, "am1, am2:", am1, am2
                print *, "basis_set_n", basis_set_n
                print *, "A, B:", A, B
                print *, "exp1, exp2", exp1, exp2
                print *, "coeff1, coeff2", coeff1, coeff2

                print *, "Got here!"
    
                call libint_overlap(am1, basis_set_n, A, exp1, coeff1, &
                                    am2, basis_set_n, B, exp2, coeff2, &
                                    overlapAB)
    
                deallocate(overlapAB)
            end do
        end do        


    end subroutine compute_S_libint



    !----------------------------------------------------------
    ! Setup the index_map for (mu,nu,lam,sig) → ERI_flat index
    !----------------------------------------------------------
    subroutine setup_eri_index_map(eri_list, n_unique, nbf, index_map)
        type(ERI_Index_type), intent(in) :: eri_list(:)
        integer*8, intent(in) :: n_unique
        integer, intent(in) :: nbf
        integer :: idx, mu, nu, lam, sig
        integer, intent(out) :: index_map(:,:,:,:)

        index_map = 0
        index_map_nbf = nbf

        do idx = 1, n_unique
            mu  = eri_list(idx)%mu
            nu  = eri_list(idx)%nu
            lam = eri_list(idx)%lam
            sig = eri_list(idx)%sig
            index_map(mu, nu, lam, sig) = idx
        end do
    end subroutine setup_eri_index_map

    !----------------------------------------------------------
    ! Main ERI computation routine using Libint
    !----------------------------------------------------------
    subroutine compute_ERI_libint(mol, shells, basis_function, n_unique, ERI_flat, eri_list, basis_set_n, index_map)
        use iso_c_binding, only : c_double, c_null_ptr
        
        type(Molecule_type),      intent(in)  :: mol
        type(Shell_type),         intent(in)  :: shells(:)
        type(BasisFunction_type), intent(in)  :: basis_function(:)
        integer*8,                intent(in)  :: n_unique
        type(ERI_Index_type),     intent(in)  :: eri_list(n_unique)
        real(c_double),           intent(out) :: ERI_flat(:)
        integer,                  intent(in)  :: basis_set_n
        integer,                  intent(out) :: index_map(:,:,:,:)
        ! Screening parameters
    

        integer                    :: s1, s2, s3, s4                     ! shell indices
        integer                    :: n1, n2, n3, n4                     ! # bf in each shell
        integer                    :: am1, am2, am3, am4, max_am
        integer                    :: i1, i2, i3, i4                     ! loop vars within shells
        integer                    :: mu, nu, lam, sig                   ! bf indices
        integer                    :: idx_val, idx_eri
        integer                    :: deriv_order, nbf
        
        type(libint_t), dimension(basis_set_n**4) :: erieval
        real(c_double), pointer    :: eri_values(:)
        real(c_double), allocatable :: F(:,:)                             ! Boys‐function buffer
        real(c_double), allocatable :: coeff1(:), coeff2(:), coeff3(:), coeff4(:)
        real(c_double), allocatable :: exp1(:),  exp2(:),  exp3(:),  exp4(:)
        real(c_double)              :: A(3), B(3), C(3), D(3)

    
        nbf = size(basis_function)
        call setup_eri_index_map(eri_list, n_unique, nbf, index_map)
    
        deriv_order = 0
        call libint2_static_init()
    
        do s1 = 1, size(shells)
            do s2 = 1, s1
                do s3 = 1, size(shells)
                    do s4 = 1, s3
    
                        n1  = size(shells(s1)%bf_idx)
                        n2  = size(shells(s2)%bf_idx)
                        n3  = size(shells(s3)%bf_idx)
                        n4  = size(shells(s4)%bf_idx)
    
                        am1 = shells(s1)%l
                        am2 = shells(s2)%l
                        am3 = shells(s3)%l
                        am4 = shells(s4)%l
                        max_am = maxval([am1, am2, am3, am4])
        
                        A      = shells(s1)%center
                        B      = shells(s2)%center
                        C      = shells(s3)%center
                        D      = shells(s4)%center
                        exp1   = shells(s1)%exponents
                        exp2   = shells(s2)%exponents
                        exp3   = shells(s3)%exponents
                        exp4   = shells(s4)%exponents
                        coeff1 = shells(s1)%coefficients
                        coeff2 = shells(s2)%coefficients
                        coeff3 = shells(s3)%coefficients
                        coeff4 = shells(s4)%coefficients
    
                        allocate(F(am1 + am2 + am3 + am4 + 1 + deriv_order, basis_set_n**4))
    
                        call build_F(am1, am2, am3, am4, deriv_order, basis_set_n, F, &
                                     [A,B,C,D], &
                                     [exp1,exp2,exp3,exp4], &
                                     [coeff1, coeff2, coeff3, coeff4])
    
                        call libint2_init_eri(erieval, max_am, c_null_ptr)
    
                        call compute_eri_f(basis_set_n, deriv_order, am1, coeff1, exp1, A, &
                                           am2, coeff2, exp2, B, &
                                           am3, coeff3, exp3, C, &
                                           am4, coeff4, exp4, D, &
                                           F, erieval)
    
                        call c_f_pointer(erieval(1)%targets(1), eri_values, [n1*n2*n3*n4])
    
                        idx_val = 0
                        do i1 = 1, n1
                            mu = shells(s1)%bf_idx(i1)
                            do i2 = 1, n2
                                nu = shells(s2)%bf_idx(i2)
                                do i3 = 1, n3
                                    lam = shells(s3)%bf_idx(i3)
                                    do i4 = 1, n4
                                        sig = shells(s4)%bf_idx(i4)
                                        idx_val = idx_val + 1
                                        idx_eri = index_map(mu, nu, lam, sig)
                                        if (idx_eri > 0) ERI_flat(idx_eri) = eri_values(idx_val)
                                    end do
                                end do
                            end do
                        end do
    
                        call libint2_cleanup_eri(erieval)
                        deallocate(F)
                    end do
                end do
            end do
        end do
    
    end subroutine compute_ERI_libint


    function compute_ERI_value_libint(mu, nu, lam, sig, shells, basis_function, basis_set_n) result(eri_val)
        use ISO_C_BINDING, ONLY: C_DOUBLE, C_NULL_PTR, C_F_POINTER
        use libint_f, ONLY: libint_t, libint2_init_eri, libint2_cleanup_eri, compute_eri_f    
        integer, intent(in) :: mu, nu, lam, sig
        type(Shell_type), intent(in) :: shells(:)
        type(BasisFunction_type), intent(in) :: basis_function(:)
        integer, intent(in) :: basis_set_n
        real(C_DOUBLE) :: eri_val
    
        integer :: mu_shell, nu_shell, lam_shell, sig_shell
        integer :: n1, n2, n3, n4, max_am, deriv_order
        integer :: am1, am2, am3, am4
        integer :: i_mu, i_nu, i_lam, i_sig
        integer :: flat_index
        TYPE(libint_t), DIMENSION(basis_set_n**4) :: erieval
        real(c_double), pointer :: eri_values(:)
        real(C_DOUBLE), allocatable :: F(:, :)
        real(c_double), allocatable :: coeff1(:), coeff2(:), coeff3(:), coeff4(:)
        real(c_double), allocatable :: exp1(:), exp2(:), exp3(:), exp4(:)
        real(c_double) :: A(3), B(3), C(3), D(3)
        real(C_DOUBLE) :: coeff_max1, coeff_max2, coeff_max3, coeff_max4, coeff_prod
        real(C_DOUBLE), parameter :: coeff_screen_tol = 1e-10  ! or as appropriate
    
        deriv_order = 0
    
        ! Get shell indices for the basis functions
        mu_shell  = basis_function(mu)%shell_idx
        nu_shell  = basis_function(nu)%shell_idx
        lam_shell = basis_function(lam)%shell_idx
        sig_shell = basis_function(sig)%shell_idx
    
        ! Get number of basis functions in each shell
        n1 = size(shells(mu_shell)%bf_idx)
        n2 = size(shells(nu_shell)%bf_idx)
        n3 = size(shells(lam_shell)%bf_idx)
        n4 = size(shells(sig_shell)%bf_idx)
    
        ! Get angular momenta
        am1 = shells(mu_shell)%l
        am2 = shells(nu_shell)%l
        am3 = shells(lam_shell)%l
        am4 = shells(sig_shell)%l
        max_am = MAXVAL([am1, am2, am3, am4])
    
        ! Get centers
        A = shells(mu_shell)%center
        B = shells(nu_shell)%center
        C = shells(lam_shell)%center
        D = shells(sig_shell)%center
    
        ! Get exponents and coefficients
        exp1 = shells(mu_shell)%exponents
        exp2 = shells(nu_shell)%exponents
        exp3 = shells(lam_shell)%exponents
        exp4 = shells(sig_shell)%exponents
    
        coeff1 = shells(mu_shell)%coefficients
        coeff2 = shells(nu_shell)%coefficients
        coeff3 = shells(lam_shell)%coefficients
        coeff4 = shells(sig_shell)%coefficients

        coeff_max1 = maxval(abs(shells(mu_shell)%coefficients))
        coeff_max2 = maxval(abs(shells(nu_shell)%coefficients))
        coeff_max3 = maxval(abs(shells(lam_shell)%coefficients))
        coeff_max4 = maxval(abs(shells(sig_shell)%coefficients))
        coeff_prod = coeff_max1 * coeff_max2 * coeff_max3 * coeff_max4
        
        if (coeff_prod < coeff_screen_tol) then
            eri_val = 0.0_C_DOUBLE
            return
        end if
        
        allocate(F(am1 + am2 + am3 + am4 + 1 + deriv_order, basis_set_n**4))
    
    
        call build_F(am1, am2, am3, am4, deriv_order, basis_set_n, F, [A,B,C,D], &
                     [exp1,exp2,exp3,exp4], [coeff1, coeff2, coeff3, coeff4])
                     
        call libint2_init_eri(erieval, max_am, C_NULL_PTR)
        call compute_eri_f(basis_set_n, deriv_order, am1, coeff1, exp1, A, &
                           am2, coeff2, exp2, B, &
                           am3, coeff3, exp3, C, &
                           am4, coeff4, exp4, D, &
                           F, erieval)
    
        eri_val = 0.0_C_DOUBLE
        if (erieval(1)%targets(1) /= C_NULL_PTR) then
            call C_F_POINTER(erieval(1)%targets(1), eri_values, [n1*n2*n3*n4])
            ! Find shell-local indices (0-based for C-order)
            i_mu  = findloc(shells(mu_shell)%bf_idx, mu, 1) - 1
            i_nu  = findloc(shells(nu_shell)%bf_idx, nu, 1) - 1
            i_lam = findloc(shells(lam_shell)%bf_idx, lam, 1) - 1
            i_sig = findloc(shells(sig_shell)%bf_idx, sig, 1) - 1
        
            if (i_mu < 0 .or. i_mu >= n1 .or. i_nu < 0 .or. i_nu >= n2 .or. &
                i_lam < 0 .or. i_lam >= n3 .or. i_sig < 0 .or. i_sig >= n4) then
                print *, "Index out of bounds!"
                print *, "mu, nu, lam, sig:", mu, nu, lam, sig
                print *, "mu_shell, nu_shell, lam_shell, sig_shell:", mu_shell, nu_shell, lam_shell, sig_shell
                print *, "i_mu, i_nu, i_lam, i_sig:", i_mu, i_nu, i_lam, i_sig
                print *, "n1, n2, n3, n4:", n1, n2, n3, n4
                stop
            end if
        
            flat_index = i_mu + n1 * (i_nu + n2 * (i_lam + n3 * i_sig)) + 1
            eri_val = eri_values(flat_index)
        end if

        if (eri_val < ERI_MIN .or. eri_val > ERI_MAX .or. abs(eri_val) < 1d-10) then
            eri_val = 0.0_C_DOUBLE
        end if
    
        
        call libint2_cleanup_eri(erieval)
        deallocate(F)
    end function compute_ERI_value_libint

    subroutine compute_ERI_libint2(mol, shells, basis_function, n_unique, ERI_flat, eri_list, basis_set_n, index_map)
        use iso_c_binding, only : c_double
        implicit none
        type(Molecule_type),      intent(in)  :: mol
        type(Shell_type),         intent(in)  :: shells(:)
        type(BasisFunction_type), intent(in)  :: basis_function(:)
        integer*8,                intent(in)  :: n_unique
        type(ERI_Index_type),     intent(in)  :: eri_list(n_unique)
        real(c_double),           intent(out) :: ERI_flat(:)
        integer,                  intent(in)  :: basis_set_n
        integer,                  intent(out) :: index_map(:,:,:,:)
    
        integer :: nbf, i, mu, nu, lam, sig
        ! Optional: for progress bar or print status
        ! integer :: report_interval
    
        ! Setup index_map (for later symmetry lookups, etc.)
        nbf = size(basis_function)
        call setup_eri_index_map(eri_list, n_unique, nbf, index_map)
    
        ! Loop over the list of unique ERI indices and compute each one
        ERI_flat = 0.0_C_DOUBLE
        do i = 1, n_unique
            mu  = eri_list(i)%mu
            nu  = eri_list(i)%nu
            lam = eri_list(i)%lam
            sig = eri_list(i)%sig
    
            ERI_flat(i) = compute_ERI_value_libint(mu, nu, lam, sig, shells, basis_function, basis_set_n)
            ! Optionally print progress
            ! if (mod(i,report_interval) == 0) print *, "ERI", i, "of", n_unique, "done"
        end do
    
    end subroutine compute_ERI_libint2


    !----------------------------------------------------------
    ! Print shell quartet integrals (for debugging)
    !----------------------------------------------------------
    SUBROUTINE print_eri(am1, am2, am3, am4, deriv_order, erieval)
        INTEGER, INTENT(IN) :: am1, am2, am3, am4, deriv_order
        TYPE(libint_t), DIMENSION(*), INTENT(IN) :: erieval
        REAL*8, DIMENSION(:), POINTER :: eri_shell_set
        INTEGER :: n1, n2, n3, n4, i_target, na, nb, nc, nd, ishell
        INTEGER, PARAMETER, DIMENSION(3) :: n_targets = [1, 12, 78]
        
        n1 = (am1 + 1)*(am1 + 2)/2
        n2 = (am2 + 1)*(am2 + 2)/2
        n3 = (am3 + 1)*(am3 + 2)/2
        n4 = (am4 + 1)*(am4 + 2)/2
        
        DO i_target = 1, n_targets(deriv_order + 1)
            WRITE (*, "(A5,1X,I2)") "Shell-set #", i_target
            CALL C_F_POINTER(erieval(1)%targets(i_target), eri_shell_set, SHAPE=[n1*n2*n3*n4])
            
            ishell = 0
            DO na = 1, n1
                DO nb = 1, n2
                    DO nc = 1, n3
                        DO nd = 1, n4
                         ishell = ishell + 1
                         WRITE (*, "(A5, I4, A11, F21.17)") &
                            "Elem ", ishell, ", (ab|cd) = ", eri_shell_set(ishell)
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE print_eri

    
end module libint_integrals