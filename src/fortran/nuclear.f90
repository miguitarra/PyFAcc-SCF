module nuclear
    use constants_struct
    use gaussian
    use overlap
    use basis_function_struct
    use molecule_struct
    use openacc
    use utils
    use io
    use boys
    implicit none

    real(8), parameter :: screen_thresholdn = 1.0D-10
    real(8), parameter :: distance_cutoffn = 10.0D0

contains

    function f(j,l,m,a,b) result(res)
        !$acc routine seq
        integer, intent(in) :: j, l, m
        real(8), intent(in) :: a, b
        integer :: k
        real(8) :: tmp, res
        res = 0.0D0
        do k = max(0,j-m), min(j,l)
            tmp = binom(l,k) * binom(m,j-k) * a**(l-k) * b**(m+k-j)
            res = res + tmp
        end do
    end function f

    function A(l,r,i,l1,l2,Ra,Rb,Rc,Rp,eps)
        !$acc routine seq
        ! -------------------------------------
        ! Factor A
        ! -------------------------------------
        !
        ! Source:
        !   Handbook of Computational Chemistry
        !   David Cook
        !   Oxford University Press
        !   1998
        !
        !--------------------------------------

        ! INPUT
        INTEGER, intent(in) :: l, r, i , l1, l2
        REAL*8, intent(in) :: Ra, Rb, Rc, Rp, eps

        ! OUTPUT
        REAL*8 :: A

        A = 1.0D0
        A = A * (-1)**(l)
        A = A * f(l,l1,l2,Rp-Ra,Rp-Rb)
        A = A * (-1)**i
        A = A * factorial(l)
        A = A * (Rp-Rc)**(l - 2*r - 2*i)
        A = A * eps**(r+i)
        A = A / factorial(r)
        A = A / factorial(i)
        A = A / factorial(l - 2*r - 2*i)

    end function A

    FUNCTION nuclear_coeff(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb,Rn,Zn) result(Vnn)
        !$acc routine seq
        ! -------------------------------------------------------------------------
        ! Compute electon-nucleus integral between two Cartesian Gaussian functions
        ! -------------------------------------------------------------------------
        !
        ! Source:
        !   Handbook of Computational Chemistry
        !   David Cook
        !   Oxford University Press
        !   1998
        !
        !--------------------------------------------------------------------------

        ! INPUT        
        INTEGER, intent(in) :: ax, ay, az, bx, by, bz   ! Angular momentum coefficients
        REAL*8, intent(in) :: aa, bb                    ! Exponential Gaussian coefficients
        REAL*8, dimension(3), intent(in) :: Ra, Rb      ! Gaussian centers
        REAL*8, dimension(3), intent(in) :: Rn          ! Nuclear position
        INTEGER, intent(in) :: Zn                       ! Nuclear charge

        ! INTERMEDIATE VARIABLES
        REAL*8 :: eps
        REAL*8 :: g                             ! Gaussian produc exponential coefficient
        REAL*8, dimension(3) :: Rp              ! Gaussian produc center
        REAL*8 :: cp                            ! Gaussian product multiplicative constant
        INTEGER :: l, r, i, m, s, j, n, t, k    ! Loop indices
        REAL*8 :: AAx, AAy, AAz                 ! Temporary calls to FUNCTION A
        INTEGER :: nu                           ! Boys function index
        real(c_double) :: vals(0:12)

        ! OUTPUT
        REAL*8 :: Vnn ! Nuclear matrix element

        CALL gaussian_product(aa,bb,Ra,Rb,g,Rp,cp)

        eps = 1.0D0 / (4.0D0 * g)

        Vnn = 0.0D0

        DO l = 0, ax+bx
            DO r = 0, FLOOR(l / 2.0)
                DO i = 0, FLOOR((l - 2*r)/2.0)
                    AAx = A(l,r,i,ax,bx,Ra(1),Rb(1),Rn(1),Rp(1),eps)

                    DO m = 0, ay+by
                        DO s = 0, FLOOR(m / 2.0)
                            DO j = 0, FLOOR((m - 2*s)/2.0)
                                AAy = A(m,s,j,ay,by,Ra(2),Rb(2),Rn(2),Rp(2),eps)

                                DO n = 0, az+bz
                                    DO t = 0, FLOOR(n / 2.0)
                                        DO k = 0, FLOOR((n - 2*t)/2.0)
                                            AAz =  A(n,t,k,az,bz,Ra(3),Rb(3),Rn(3),Rp(3),eps)

                                            nu = l + m + n - 2 * (r + s + t) - (i + j + k)

                                            Vnn = Vnn + AAx * AAy * AAz * boys(nu,g*DOT_PRODUCT(Rp-Rn,Rp-Rn))
                                        END DO ! k
                                    END DO ! t
                                END DO ! n
                            END DO ! j
                        END DO ! s
                    END DO ! m
                END DO ! i
            END DO ! r
        END DO ! l

        Vnn = Vnn * (-Zn) * norm(ax,ay,az,aa) * norm(bx,by,bz,bb) * cp * 2.0D0 * PI / g

    END FUNCTION nuclear_coeff

    subroutine V_nuclear(Kf, c, basis_functions, Rn, Zn, Nn, V)
        use iso_fortran_env
        implicit none
    
        integer, intent(in) :: Kf, c, Nn
        type(BasisFunction_type), intent(in) :: basis_functions(Kf)
        real(8), intent(in) :: Rn(3, Nn)
        integer, intent(in) :: Zn(Nn)
        real(8), intent(out) :: V(Kf, Kf)
    
        ! Local copies of basis information
        real(8) :: basis_D(Kf, c), basis_A(Kf, c)
        integer :: basis_L(Kf, 3)
        real(8) :: basis_R(Kf, 3)
    
        ! Loop variables
        integer :: i, j, k, l, n
        real(8) :: dprod, dist2, tmp
        real(8), dimension(3) :: Ri, Rj
    
        ! 3D accumulation array
        real(8), allocatable :: Vtmp(:,:,:)
        real(8) :: V_local(Kf, Kf)
    
        ! Initialize output
        V = 0.0d0
        V_local = 0.0d0
        allocate(Vtmp(Nn, Kf, Kf))
        Vtmp = 0.0d0
    
        ! Extract basis data into plain arrays
        do i = 1, Kf
            basis_D(i,:) = basis_functions(i)%coefficients
            basis_A(i,:) = basis_functions(i)%exponents
            basis_L(i,:) = basis_functions(i)%ang_mom
            basis_R(i,:) = basis_functions(i)%center
        end do
    
        !$acc data copyin(basis_D, basis_A, basis_L, basis_R, Rn, Zn), create(Vtmp) async(3)
        !$acc parallel loop collapse(3) private(Ri, Rj, k, l, dprod, dist2, tmp) default(present) async(3)
        do n = 1, Nn
            do i = 1, Kf
                do j = 1, Kf
                    Ri = basis_R(i,:)
                    Rj = basis_R(j,:)
                    dist2 = sum((Ri - Rj)**2)
    
                    if (dist2 > distance_cutoffn**2) cycle
    
                    tmp = 0.0d0
                    do k = 1, c
                        do l = 1, c
                            dprod = basis_D(i,k) * basis_D(j,l)
                            if (abs(dprod) < screen_thresholdn) cycle
    
                            tmp = tmp + dprod * nuclear_coeff( &
                                basis_L(i,1), basis_L(i,2), basis_L(i,3), &
                                basis_L(j,1), basis_L(j,2), basis_L(j,3), &
                                basis_A(i,k), basis_A(j,l), &
                                Ri, Rj, Rn(:,n), Zn(n) )
                        end do
                    end do
                    Vtmp(n, i, j) = tmp
                end do
            end do
        end do
        !$acc end parallel loop
    
        !$acc update self(Vtmp) async(3)
        !$acc wait(3)
    
        ! Reduce Vtmp over n to get V_local
        V_local = 0.0d0
        do n = 1, Nn
            V_local = V_local + Vtmp(n, :, :)
        end do
    
        V = V + V_local
        !$acc end data
    
        deallocate(Vtmp)
    
    end subroutine V_nuclear

    
end module nuclear
