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
    use cartesian
    use basis_struct
    use openacc

    implicit none
    
    real(8), parameter :: screen_threshold = 1.0D-10
    real(8), parameter :: distance_cutoff = 10.0D0

contains

    function Si(a,b,aa,bb,Rai,Rbi,Ri)
        !$acc routine seq
        ! ------------------------------------------------------------------------------------------------
        ! Compute overlap integral between two unnormalized Cartesian Gaussian functions along direction i
        ! ------------------------------------------------------------------------------------------------

        !$acc routine seq

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
        !$acc routine seq
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

        !$acc routine seq

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
        !$acc routine seq
        real(8), intent(in) :: aa, bb, Rai, Rbi
        real(8) :: S00
        S00 = DSQRT(PI / (aa + bb)) * DEXP(- (aa * bb) / (aa + bb) * (Rai - Rbi)**2.0D0)
    end function S00

    !-----------------------------------------------------------
    ! OS recursion for overlap integrals
    !-----------------------------------------------------------
    function Sij(a, b, aa, bb, Rai, Rbi, Rpi)
        !$acc routine seq
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
        !$acc routine seq
    
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

    subroutine compute_overlap_shell_pair(xA, xB, yA, yB, zA, zB, lA, lB, startA, endA, startB, endB, exponents, coefficients, S_block)
        !$acc routine seq
        ! Inputs
        integer, intent(in) :: lA, lB, startA, endA, startB, endB
        real(wp), intent(in) :: xA, xB, yA, yB, zA, zB
        real(wp), intent(in) :: exponents(:), coefficients(:)
    
        ! Output
        real(wp), intent(out) :: S_block(:,:)   ! (nA, nB)
    
        ! Local variables
        integer :: nA, nB
        integer :: muA, muB
        integer :: axA, ayA, azA, axB, ayB, azB
        integer :: kA, kB
        real(wp) :: tmp
        real(wp), dimension(3) :: RA, RB
    
        ! Fixed-size arrays for angular momenta
        integer, parameter :: l_max = 3               ! adjust to max l in your shells
        integer, parameter :: nAmax = (l_max+1)*(l_max+2)/2
        integer, parameter :: nBmax = (l_max+1)*(l_max+2)/2
        integer :: angmomA(nAmax,3)
        integer :: angmomB(nBmax,3)
    
        ! Number of basis functions in each shell
        nA = (lA + 1)*(lA + 2)/2
        nB = (lB + 1)*(lB + 2)/2
    
        ! Shell centers
        RA = [xA, yA, zA]
        RB = [xB, yB, zB]
    
        ! Generate Cartesian angular momenta
        call generate_cartesian_angmoms(lA, angmomA)
        call generate_cartesian_angmoms(lB, angmomB)
    
        ! Loop over basis function pairs
        do muA = 1, nA
            axA = angmomA(muA,1)
            ayA = angmomA(muA,2)
            azA = angmomA(muA,3)
            do muB = 1, nB
                axB = angmomB(muB,1)
                ayB = angmomB(muB,2)
                azB = angmomB(muB,3)
                tmp = 0.0_wp
                do kA = startA, endA
                    do kB = startB, endB
                        tmp = tmp + coefficients(kA) * coefficients(kB) &
                              * overlap_coeff_OS(axA, ayA, azA, axB, ayB, azB, &
                                                 exponents(kA), exponents(kB), RA, RB)
                    end do
                end do
                S_block(muA, muB) = tmp
            end do
        end do
    
    end subroutine compute_overlap_shell_pair


    !-----------------------------------------------------------
    ! Compute the overall overlap matrix with screening and distance cutoff.
    !-----------------------------------------------------------
    subroutine S_overlap(mol, shells, S)
        ! input
        type(Molecule_type), intent(in) :: mol
        type(Shell_type), intent(in) :: shells(:)
        real(wp), intent(out) :: S(:,:)
        
        ! local
        integer :: n_shells, sA, sB, nA, nB, muA, muB, gA, gB, idx, i, j
        integer :: prim_index, coeffA_start, coeffA_end, coeffB_start, coeffB_end
        integer, allocatable :: shell_start(:)
        integer, allocatable :: prim_start(:)
        integer :: total_prim
        real(wp), allocatable :: S_block(:,:)
        real(wp), allocatable :: shell_x(:), shell_y(:), shell_z(:), shell_coeff(:), shell_exp(:)
        integer, allocatable :: shell_prim(:), shell_l(:)
        integer :: n_upper
    
        n_shells = size(shells)
        allocate(shell_start(n_shells+1))
        allocate(shell_prim(n_shells), shell_l(n_shells), shell_x(n_shells), shell_y(n_shells), shell_z(n_shells))
    
        ! compute total number of primitives and prim_start array
        total_prim = 0
        do sA = 1, n_shells
            total_prim = total_prim + shells(sA)%num_prim
        end do
        allocate(shell_coeff(total_prim), shell_exp(total_prim))
        allocate(prim_start(n_shells+1))
    
        ! fill arrays to make it OpenACC compatible
        shell_start(1) = 1
        prim_index = 1
        prim_start(1) = 1
        do sA = 1, n_shells
            shell_x(sA) = shells(sA)%x
            shell_y(sA) = shells(sA)%y
            shell_z(sA) = shells(sA)%z
            shell_prim(sA) = shells(sA)%num_prim
            shell_l(sA) = shells(sA)%l_num
    
            ! correct slice: prim_index : prim_index + num_prim - 1
            shell_coeff(prim_index : prim_index + shell_prim(sA) - 1) = shells(sA)%coefficients
            shell_exp  (prim_index : prim_index + shell_prim(sA) - 1) = shells(sA)%exponents
    
            prim_index = prim_index + shell_prim(sA)
            prim_start(sA+1) = prim_index
    
            nA = (shell_l(sA) + 1)*(shell_l(sA) + 2)/2
            shell_start(sA+1) = shell_start(sA) + nA
        end do
    
        S(:,:) = 0.0_wp
        n_upper = n_shells * (n_shells + 1) / 2
    
        
        !$acc data copyin(shell_x, shell_y, shell_z, shell_l, shell_coeff, shell_prim, shell_exp, shell_start, prim_start)
        !$acc parallel loop private(sA, sB, nA, nB, coeffA_start, coeffA_end, coeffB_start, coeffB_end, muA, muB, gA, gB) default(present)
        do idx = 1, n_upper
            call idx_upper_triangular(n_shells, idx, sA, sB)
            nA = (shell_l(sA) + 1)*(shell_l(sA) + 2)/2
            nB = (shell_l(sB) + 1)*(shell_l(sB) + 2)/2
    
            coeffA_start = prim_start(sA)
            coeffA_end   = prim_start(sA+1) - 1
            coeffB_start = prim_start(sB)
            coeffB_end   = prim_start(sB+1) - 1
    
            allocate(S_block(nA, nB))
    
            call compute_overlap_shell_pair(shell_x(sA), shell_x(sB), &
                                           shell_y(sA), shell_y(sB), &
                                           shell_z(sA), shell_z(sB), &
                                           shell_l(sA), shell_l(sB), &
                                           coeffA_start, coeffA_end, &
                                           coeffB_start, coeffB_end, &
                                           shell_exp, shell_coeff, S_block)
    
            do muA = 1, nA
                gA = shell_start(sA) + muA - 1
                do muB = 1, nB
                    gB = shell_start(sB) + muB - 1
    
                    S(gA, gB) = S_block(muA, muB)
                end do
            end do
    
            deallocate(S_block)
    
        end do
        !$acc end parallel loop
        !$acc update host(S)
    
        do i = 1, size(S, 1)
            do j = i+1, size(S, 1)
                S(j, i) = S(i, j)
            end do
        end do

        !$acc end data
        
    end subroutine S_overlap




end module overlap
