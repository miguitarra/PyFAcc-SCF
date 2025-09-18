! --------------------------------------------------------------------
! Taken from github project: https://github.com/RMeli/Hartree-Fock/tree/master
! ---------------------------------------------------------------------
module kinetic
    use constants_struct
    use utils
    use gaussian
    use overlap
    use basis_struct
    use openacc

    implicit none

    ! Thresholds
    real(8), parameter :: screen_thresholdk = 1.0D-10
    real(8), parameter :: distance_cutoffk = 10.0D0

contains

    function Ki(ac,a1,a2,bc,b1,b2,aa,bb,Ra,Rb,Ra1,Rb1,Ra2,Rb2,Rc,R1,R2,cp)
        !$acc routine seq
        integer, intent(in) :: ac, a1, a2, bc, b1, b2
        real(8), intent(in) :: aa, bb
        real(8), intent(in) :: Ra, Rb, Ra1, Rb1, Ra2, Rb2, Rc, R1, R2, cp
        real(8) :: k, Ki

        k = 0.0D0
        k = k + ac * bc * Si(ac - 1, bc - 1, aa, bb, Ra, Rb, Rc)
        k = k - 2.0D0 * aa * bc * Si(ac + 1, bc - 1, aa, bb, Ra, Rb, Rc)
        k = k - 2.0D0 * ac * bb * Si(ac - 1, bc + 1, aa, bb, Ra, Rb, Rc)
        k = k + 4.0D0 * aa * bb * Si(ac + 1, bc + 1, aa, bb, Ra, Rb, Rc)
        k = 0.5D0 * k

        Ki = cp * (PI / (aa + bb))**1.5D0 * k
        Ki = Ki * Si(a1, b1, aa, bb, Ra1, Rb1, R1)
        Ki = Ki * Si(a2, b2, aa, bb, Ra2, Rb2, R2)
    end function Ki

    function kinetic_coeff(ax, ay, az, bx, by, bz, aa, bb, Ra, Rb) result(K)
        !$acc routine seq
        integer, intent(in) :: ax, ay, az, bx, by, bz
        real(8), intent(in) :: aa, bb
        real(8), dimension(3), intent(in) :: Ra, Rb
        real(8) :: pp, cp, K
        real(8), dimension(3) :: Rp

        call gaussian_product(aa, bb, Ra, Rb, pp, Rp, cp)

        K = 0.0D0
        K = K + Ki(ax, ay, az, bx, by, bz, aa, bb, Ra(1), Rb(1), Ra(2), Rb(2), Ra(3), Rb(3), Rp(1), Rp(2), Rp(3), cp)
        K = K + Ki(ay, az, ax, by, bz, bx, aa, bb, Ra(2), Rb(2), Ra(3), Rb(3), Ra(1), Rb(1), Rp(2), Rp(3), Rp(1), cp)
        K = K + Ki(az, ax, ay, bz, bx, by, aa, bb, Ra(3), Rb(3), Ra(1), Rb(1), Ra(2), Rb(2), Rp(3), Rp(1), Rp(2), cp)
        K = K * norm(ax, ay, az, aa) * norm(bx, by, bz, bb)
    end function kinetic_coeff

    subroutine compute_kinetic_shell_pair(xA, xB, yA, yB, zA, zB, lA, lB, startA, endA, startB, endB, exponents, coefficients, T_block)
        !$acc routine seq
        ! Inputs
        integer, intent(in) :: lA, lB, startA, endA, startB, endB
        real(wp), intent(in) :: xA, xB, yA, yB, zA, zB
        real(wp), intent(in) :: exponents(:), coefficients(:)
    
        ! Output
        real(wp), intent(out) :: T_block(:,:)   ! (nA, nB)
    
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
                              * kinetic_coeff(axA, ayA, azA, axB, ayB, azB, &
                                                 exponents(kA), exponents(kB), RA, RB)
                    end do
                end do
                T_block(muA, muB) = tmp
            end do
        end do
    
    end subroutine compute_kinetic_shell_pair
    

    subroutine T_kinetic(mol, shells, T)
        type(Molecule_type), intent(in) :: mol
        type(Shell_type), intent(in) :: shells(:)
        real(8), intent(out) :: T(:, :)
    
        integer :: n_shells, sA, sB, nA, nB, muA, muB, gA, gB, idx, i, j
        integer :: prim_index, coeffA_start, coeffA_end, coeffB_start, coeffB_end
        integer, allocatable :: shell_start(:)
        integer, allocatable :: prim_start(:)
        integer :: total_prim
        real(wp), allocatable :: T_block(:,:)
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
    
        T(:,:) = 0.0_wp
        n_upper = n_shells * (n_shells + 1) / 2
    
        
        !$acc data copyin(shell_x, shell_y, shell_z, shell_l, shell_coeff, shell_prim, shell_exp, shell_start, prim_start) async(2)
        !$acc parallel loop private(sA, sB, nA, nB, coeffA_start, coeffA_end, coeffB_start, coeffB_end, muA, muB, gA, gB) default(present) async(2)
        do idx = 1, n_upper
            call idx_upper_triangular(n_shells, idx, sA, sB)
            nA = (shell_l(sA) + 1)*(shell_l(sA) + 2)/2
            nB = (shell_l(sB) + 1)*(shell_l(sB) + 2)/2
    
            coeffA_start = prim_start(sA)
            coeffA_end   = prim_start(sA+1) - 1
            coeffB_start = prim_start(sB)
            coeffB_end   = prim_start(sB+1) - 1
    
            allocate(T_block(nA, nB))
    
            call compute_kinetic_shell_pair(shell_x(sA), shell_x(sB), &
                                           shell_y(sA), shell_y(sB), &
                                           shell_z(sA), shell_z(sB), &
                                           shell_l(sA), shell_l(sB), &
                                           coeffA_start, coeffA_end, &
                                           coeffB_start, coeffB_end, &
                                           shell_exp, shell_coeff, T_block)
    
            do muA = 1, nA
                gA = shell_start(sA) + muA - 1
                do muB = 1, nB
                    gB = shell_start(sB) + muB - 1
    
                    T(gA, gB) = T_block(muA, muB)
                end do
            end do
    
            deallocate(T_block)
    
        end do
        !$acc end parallel loop 
        !$acc update host(T)
    
        do i = 1, size(T, 1)
            do j = i+1, size(T, 1)
                T(j, i) = T(i, j)
            end do
        end do

        !$acc end data
    end subroutine T_kinetic


end module kinetic
