! --------------------------------------------------------------------
!
! Taken from GH project https://github.com/RMeli/Hartree-Fock
!
! ---------------------------------------------------------------------

module overlap

    use constants_struct
    use utils
    use gaussian
    use la
    use basis_function_struct
    use openacc

    implicit none
    
    real(8), parameter :: screen_threshold = 1.0D-12
    real(8), parameter :: distance_cutoff = 10.0D0

contains

    function Si(a,b,aa,bb,Rai,Rbi,Ri)
        ! ------------------------------------------------------------------------------------------------
        ! Compute overlap integral between two unnormalized Cartesian Gaussian functions along direction i
        ! ------------------------------------------------------------------------------------------------

        !acc routine seq

        ! INPUT
        integer, intent(in) :: a,b ! Angular momentum coefficients of the Gaussians along direction i
        real*8, intent(in) :: aa, bb ! Exponential coefficients of the Gaussians
        real*8, intent(in) :: Rai, Rbi ! Centers of the Gaussians
        real*8, intent(in) :: Ri ! Center of the Gaussians product

        ! INTERMEDIATE VARIABLE
        real*8 :: tmp
        integer :: i ! Loop index
        integer :: j ! Loop index

        ! OUTPUT
        real*8 :: Si

        Si = 0.0D0

        do i = 0, a
            do j = 0, b
                IF (MOD(i+j,2) == 0) THEN
                    tmp = binom(a,i) * binom(b,j) * factorial2(i + j - 1)
                    tmp = tmp * (Ri-Rai)**(a-i)
                    tmp = tmp * (Ri-Rbi)**(b-j)
                    tmp = tmp / (2.0D0 * (aa + bb))**((i + j) / 2.0D0)

                    Si = Si + tmp
                end IF
            end do
        end do
    end function Si



    function overlap_coeff(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb) result(S)
        ! -----------------------------------------------------------------
        ! Compute overlap integral between two Cartesian Gaussian functions
        ! -----------------------------------------------------------------
        !
        ! Source:
        !   The Mathematica Journal
        !   Evaluation of Gaussian Molecular Integrals
        !   I. Overlap Integrals
        !   Minhhuy Hô and Julio Manuel Hernández-Pérez
        !
        ! -----------------------------------------------------------------

        !acc routine seq

        ! INPUT
        integer, intent(in) :: ax, ay, az, bx, by, bz ! Angular momentum coefficients
        real*8, intent(in) :: aa, bb ! Exponential Gaussian coefficients
        real*8, dimension(3), intent(in) :: Ra, Rb ! Gaussian centers

        ! INTERMEDIATE VARIABLES
        real*8 :: pp ! Gaussian produc exponential coefficient
        real*8, dimension(3) :: Rp ! Gaussian produc center
        real*8 :: cp ! Gaussian product multiplicative constant

        ! OUTPUT
        real*8 :: S

        CALL gaussian_product(aa,bb,Ra,Rb,pp,Rp,cp) ! Compute PP, RP and CP

        S = 1
        S = S * Si(ax,bx,aa,bb,Ra(1),Rb(1),Rp(1))       ! Overlap along x
        S = S * Si(ay,by,aa,bb,Ra(2),Rb(2),Rp(2))       ! Overlap along y
        S = S * Si(az,bz,aa,bb,Ra(3),Rb(3),Rp(3))       ! Overlap along z
        S = S * norm(ax,ay,az,aa) * norm(bx,by,bz,bb)   ! Normalization of Gaussian functions
        S = S * cp                                      ! Gaussian product factor
        S = S * (PI / pp)**(3./2.)                      ! Solution of Gaussian integral

    end function overlap_coeff

    ! ------------------
    ! OBARA-SAIKA SCHEME
    ! ------------------
    !-----------------------------------------------------------
    ! Overlap for s-type Gaussians (zero angular momentum)
    !-----------------------------------------------------------
    function S00(aa, bb, Rai, Rbi)
        real(8), intent(in) :: aa, bb, Rai, Rbi
        real(8) :: S00
        S00 = DSQRT(PI / (aa + bb)) * DEXP(- (aa * bb) / (aa + bb) * (Rai - Rbi)**2.0D0)
    end function S00

    !-----------------------------------------------------------
    ! OS recursion for overlap integrals
    !-----------------------------------------------------------
    function Sij(a, b, aa, bb, Rai, Rbi, Rpi)
        integer, intent(in) :: a, b
        real(8), intent(in) :: aa, bb, Rai, Rbi, Rpi
        real(8), dimension(-1:a+1, -1:b+1) :: S
        integer :: i, j
        real(8) :: Sij

        S(-1, :) = 0.0D0
        S(:, -1) = 0.0D0
        S(0, 0) = S00(aa, bb, Rai, Rbi)

        do i = 0, a
            do j = 0, b
                if (i + 1 <= a) then
                    S(i+1, j) = (Rpi - Rai) * S(i, j) + 1.0D0 / (2.0D0 * (aa + bb)) * (i * S(i-1, j) + j * S(i, j-1))
                end if
                if (j + 1 <= b) then
                    S(i, j+1) = (Rpi - Rbi) * S(i, j) + 1.0D0 / (2.0D0 * (aa + bb)) * (i * S(i-1, j) + j * S(i, j-1))
                end if
            end do
        end do

        Sij = S(a, b)
    end function Sij

    !-----------------------------------------------------------
    ! Fast overlap coefficient using OS recursion along each direction.
    !-----------------------------------------------------------
    function overlap_coeff_OS(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb) result(S)
        ! ------------------------------------------------------------------------------------
        ! Compute overlap integral between two Cartesian Gaussian functions using OS recursion
        ! ------------------------------------------------------------------------------------

        !acc routine seq

        ! INPUT
        INTEGER, intent(in) :: ax, ay, az, bx, by, bz ! Angular momentum coefficients
        REAL*8, intent(in) :: aa, bb ! Exponential Gaussian coefficients
        REAL*8, dimension(3), intent(in) :: Ra, Rb ! Gaussian centers

        ! INTERMEDIATE VARIABLES
        REAL*8 :: pp ! Gaussian produc exponential coefficient
        REAL*8, dimension(3) :: Rp ! Gaussian produc center
        REAL*8 :: cp ! Gaussian product multiplicative constant

        ! OUTPUT
        REAL*8 :: S

        CALL gaussian_product(aa,bb,Ra,Rb,pp,Rp,cp) ! Compute PP, RP and CP

        S = 1
        S = S * Sij(ax,bx,aa,bb,Ra(1),Rb(1),Rp(1))       ! Overlap along x
        S = S * Sij(ay,by,aa,bb,Ra(2),Rb(2),Rp(2))       ! Overlap along y
        S = S * Sij(az,bz,aa,bb,Ra(3),Rb(3),Rp(3))       ! Overlap along z
        S = S * norm(ax,ay,az,aa) * norm(bx,by,bz,bb)    ! Normalization of Gaussian functions

    end function overlap_coeff_OS

    !-----------------------------------------------------------
    ! Compute the overall overlap matrix with screening and distance cutoff.
    !-----------------------------------------------------------
    subroutine S_overlap(Kf, c, basis_functions, S)
    
        use openacc
        implicit none
    
        integer, intent(in) :: Kf, c
        type(BasisFunction_type), intent(in) :: basis_functions(Kf)
        real(8), intent(out) :: S(Kf, Kf)
    
        real(8) :: basis_D(Kf, c), basis_A(Kf, c)
        integer :: basis_L(Kf, 3)
        real(8) :: basis_R(Kf, 3)
    
        integer :: i, j, k, l
        real(8) :: tmp, dprod, dist2
        real(8), dimension(3) :: Ri, Rj
        integer :: idx, n_upper
        real(8) :: f
    
        ! Precompute number of upper-triangle elements
        n_upper = Kf * (Kf + 1) / 2
    
        ! Copy basis data into plain arrays (OpenACC-friendly)
        do i = 1, Kf
            basis_D(i, :) = basis_functions(i)%coefficients
            basis_A(i, :) = basis_functions(i)%exponents
            basis_L(i, :) = basis_functions(i)%ang_mom
            basis_R(i, :) = basis_functions(i)%center
        end do
    
        S = 0.0D0
    
        !$acc data copyin(basis_D, basis_A, basis_L, basis_R) copyout(S)
        !$acc parallel loop gang vector private(idx, i, j, k, l, tmp, dprod, dist2, Ri, Rj, f)
        do idx = 1, n_upper
    
            ! Map flat index idx to (i, j) in upper triangle (1-based)
            f = 0.5d0 * (-1.0d0 + sqrt(1.0d0 + 8.0d0 * dble(idx)))
            i = int(f)
            if (f /= dble(i)) i = i + 1
            j = idx - (i - 1) * i / 2
    
            Ri = basis_R(i, :)
            Rj = basis_R(j, :)
            dist2 = (Ri(1) - Rj(1))**2 + (Ri(2) - Rj(2))**2 + (Ri(3) - Rj(3))**2
    
            if (dist2 > distance_cutoff**2) cycle
    
            tmp = 0.0D0
            do k = 1, c
                do l = 1, c
                    dprod = basis_D(i, k) * basis_D(j, l)
                    if (abs(dprod) < screen_threshold) cycle
    
                    tmp = tmp + dprod * overlap_coeff_OS( &
                        basis_L(i,1), basis_L(i,2), basis_L(i,3), &
                        basis_L(j,1), basis_L(j,2), basis_L(j,3), &
                        basis_A(i, k), basis_A(j, l), &
                        Ri, Rj )
                end do
            end do
    
            S(i, j) = tmp
            if (i /= j) S(j, i) = tmp
        end do
        !$acc end parallel loop
        !$acc end data
    
    end subroutine S_overlap


end module overlap
