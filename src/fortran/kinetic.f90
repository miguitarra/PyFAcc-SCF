! --------------------------------------------------------------------
! Taken from github project: https://github.com/RMeli/Hartree-Fock/tree/master
! ---------------------------------------------------------------------
module kinetic
    use constants_struct
    use utils
    use gaussian
    use overlap
    use basis_function_struct
    use molecule_struct
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

    subroutine T_kinetic(Kf, c, basis_functions, T)
        integer, intent(in) :: Kf, c
        type(BasisFunction_type), intent(in) :: basis_functions(Kf)
        real(8), intent(out) :: T(Kf, Kf)
    
        real(8) :: basis_D(Kf, c), basis_A(Kf, c)
        integer :: basis_L(Kf, 3)
        real(8) :: basis_R(Kf, 3)
    
        integer :: i, j, k, l
        real(8) :: tmp, dprod, dist2
        real(8) :: Ri(3), Rj(3)
    
        ! Extract basis data
        do i = 1, Kf
            basis_D(i,:) = basis_functions(i)%coefficients
            basis_A(i,:) = basis_functions(i)%exponents
            basis_L(i,:) = basis_functions(i)%ang_mom
            basis_R(i,:) = basis_functions(i)%center
        end do
        
        !$acc data copyin(basis_D, basis_A, basis_L, basis_R) async(2)
        !$acc parallel loop collapse(2) private(k, l, tmp, dprod, Ri, Rj, dist2) default(present) async(2)
        do i = 1, Kf
            do j = 1, Kf
                Ri = basis_R(i, :)
                Rj = basis_R(j, :)
                dist2 = sum((Ri(:) - Rj(:))**2)
    
                if (dist2 > distance_cutoffk**2) cycle
    
                tmp = 0.0D0
                do k = 1, c
                    do l = 1, c
                        dprod = basis_D(i,k) * basis_D(j,l)
                        if (abs(dprod) < screen_thresholdk) cycle
    
                        tmp = tmp + dprod * kinetic_coeff( &
                              basis_L(i,1), basis_L(i,2), basis_L(i,3), &
                              basis_L(j,1), basis_L(j,2), basis_L(j,3), &
                              basis_A(i,k), basis_A(j,l), &
                              Ri, Rj )
                    end do
                end do
                T(i,j) = tmp
            end do
        end do
        !$acc end parallel loop
        !$acc update self(T) async(2)
        !$acc end data
    end subroutine T_kinetic


end module kinetic
