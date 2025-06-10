module integrals
    use ISO_C_BINDING, ONLY: C_DOUBLE, C_F_POINTER, C_NULL_PTR
    USE libint_f, ONLY: libint_t, libint2_static_init, libint2_static_cleanup, libint2_build, libint2_max_am_eri, &
                        compute_eri_f, libint2_init_eri, libint2_cleanup_eri

    use molecule_struct
    use basis_function_struct
    use electronic
    use utils
    use boys

    implicit none

    integer, allocatable, public :: index_map(:,:,:,:)
    integer, private :: index_map_nbf = 0

contains

    !----------------------------------------------------------
    ! Setup the index_map for (mu,nu,lam,sig) → ERI_flat index
    !----------------------------------------------------------
    subroutine setup_eri_index_map(eri_list, n_unique, nbf)
        type(ERI_Index_type), intent(in) :: eri_list(:)
        integer*8, intent(in) :: n_unique
        integer, intent(in) :: nbf
        integer :: idx, mu, nu, lam, sig

        if (allocated(index_map)) deallocate(index_map)
        allocate(index_map(nbf, nbf, nbf, nbf))
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
    subroutine compute_ERI_libint(mol, shells, basis_function, n_unique, ERI_flat, eri_list, basis_set_n)
        type(Molecule_type), intent(in) :: mol
        type(Shell_type), intent(in) :: shells(:)
        type(BasisFunction_type), intent(in) :: basis_function(:)
        integer*8, intent(in) :: n_unique
        type(ERI_Index_type), intent(in) :: eri_list(n_unique)
        real(c_double), intent(out) :: ERI_flat(:)
        integer, intent(in) :: basis_set_n

        ! Locals
        integer :: s1, s2, s3, s4, n1, n2, n3, n4, max_am, i, s, deriv_order, nbf
        integer :: i1, i2, i3, i4, mu, nu, lam, sig, idx_val, idx_eri
        integer :: am1, am2, am3, am4
        TYPE(libint_t), DIMENSION(basis_set_n**4) :: erieval
        real(c_double), pointer :: eri_values(:)
        real(C_DOUBLE), allocatable :: F(:, :)
        real(c_double), allocatable :: coeff1(:), coeff2(:), coeff3(:), coeff4(:)
        real(c_double), allocatable :: exp1(:), exp2(:), exp3(:), exp4(:)
        real(c_double) :: A(3), B(3), C(3), D(3)

        nbf = size(basis_function)

        ! Ensure index_map is set up for this nbf!
        if (.not.allocated(index_map) .or. index_map_nbf /= nbf) then
            call setup_eri_index_map(eri_list, n_unique, nbf)
        end if

        deriv_order = 0

        do s1 = 1, size(shells)
            do s2 = 1, s1
                do s3 = 1, size(shells)
                    do s4 = 1, s3

                        ! Get shell data
                        n1 = size(shells(s1)%bf_idx)
                        n2 = size(shells(s2)%bf_idx)
                        n3 = size(shells(s3)%bf_idx)
                        n4 = size(shells(s4)%bf_idx)

                        am1 = shells(s1)%l
                        am2 = shells(s2)%l
                        am3 = shells(s3)%l
                        am4 = shells(s4)%l
                        max_am = MAXVAL([am1, am2, am3, am4])

                        ! Prepare input for Libint
                        A = shells(s1)%center
                        B = shells(s2)%center
                        C = shells(s3)%center
                        D = shells(s4)%center
                        exp1 = shells(s1)%exponents
                        exp2 = shells(s2)%exponents
                        exp3 = shells(s3)%exponents
                        exp4 = shells(s4)%exponents
                        coeff1 = shells(s1)%coefficients
                        coeff2 = shells(s2)%coefficients
                        coeff3 = shells(s3)%coefficients
                        coeff4 = shells(s4)%coefficients

                        allocate(F(am1 + am2 + am3 + am4 + 1 + deriv_order, basis_set_n**4))

                        call build_F(am1, am2, am3, am4, deriv_order, basis_set_n, F, [A,B,C,D], &
                                     [exp1,exp2,exp3,exp4], [coeff1, coeff2, coeff3, coeff4])

                        call libint2_static_init()
                        call libint2_init_eri(erieval, max_am, C_NULL_PTR)

                        call compute_eri_f(basis_set_n, deriv_order, am1, coeff1, exp1, A, &
                                           am2, coeff2, exp2, B, &
                                           am3, coeff3, exp3, C, &
                                           am4, coeff4, exp4, D, &
                                           F, erieval)
                        call C_F_POINTER(erieval(1)%targets(1), eri_values, [n1*n2*n3*n4])

                        ! Unpack only needed integrals using index_map
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
                                        ! Use index_map for fast lookup
                                        idx_eri = index_map(mu, nu, lam, sig)
                                        if (idx_eri > 0) then
                                            ERI_flat(idx_eri) = eri_values(idx_val)
                                        end if
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

    subroutine libint_calc_1body_ints(n_shells, n_basis, centers, exponents, coefficients, l_shells, shell_offsets, n_primitives, result, op_type) bind(C)
        integer(c_int), value :: n_shells, n_basis, op_type
        real(c_double), intent(in) :: centers(3*n_shells)
        real(c_double), intent(in) :: exponents(*)
        real(c_double), intent(in) :: coefficients(*)
        integer(c_int), intent(in) :: l_shells(n_shells)
        integer(c_int), intent(in) :: shell_offsets(n_shells) ! -> offsets of the flattened array 
        integer(c_int), intent(in) :: n_primitives(n_shells)
        real(c_double), intent(out) :: result(n_basis*n_basis)
    end subroutine

end module integrals