module boys
    use iso_c_binding
    use utils
    implicit none

contains


    pure FUNCTION boys(nu, x)
        !$acc routine seq
        INTEGER, INTENT(IN) :: nu
        REAL*8, INTENT(IN)  :: x
        REAL*8              :: boys
    
        ! Local variables
        INTEGER :: i, N
        REAL*8 :: t, dt, sum, exponent
    
        boys = 0.0D0
    
        IF (x .LE. 1.0D-6) THEN
            ! Taylor expansion for small x
            boys = 1.0D0 / (2.0D0 * nu + 1.0D0) - x / (2.0D0 * nu + 3.0D0)
        ELSE
            ! Numerical integration using trapezoidal rule
            N = 1000
            dt = 1.0D0 / N
            sum = 0.0D0
    
            DO i = 0, N
                t = i * dt
                exponent = EXP(-x * t * t)
                IF (i == 0 .OR. i == N) THEN
                    sum = sum + 0.5D0 * t**(2 * nu) * exponent
                ELSE
                    sum = sum + t**(2 * nu) * exponent
                END IF
            END DO
    
            boys = sum * dt
        END IF
    
    END FUNCTION boys


    subroutine build_F(am1, am2, am3, am4, deriv_order, contrdepth, F, centers, alphas, coeffs)
        !$acc routine seq
        implicit none
        integer, intent(in) :: am1, am2, am3, am4, deriv_order, contrdepth
        real(C_DOUBLE), intent(out) :: F(:,:)
        real(8), intent(in) :: centers(4,3), alphas(4,contrdepth), coeffs(4, contrdepth)
    
        integer :: nF, Np
        integer :: p1, p2, p3, p4, p1234
        integer, value :: n
        real(C_DOUBLE), value :: x

    
        nF = am1 + am2 + am3 + am4 + 1 + deriv_order
        Np = contrdepth**4
    
        F = 0.0d0
    
        do p1 = 1, contrdepth
        do p2 = 1, contrdepth
          do p3 = 1, contrdepth
            do p4 = 1, contrdepth
              p1234 = (p1-1)*contrdepth**3 + (p2-1)*contrdepth**2 + (p3-1)*contrdepth + p4
    
              ! Compute Boys function argument x
              x = compute_boys_x(centers(1,:), centers(2,:), centers(3,:), centers(4,:), &
                                 alphas(1,:), coeffs(1,:), alphas(2,:), coeffs(2,:), &
                                 alphas(3,:), coeffs(3,:), alphas(4,:), coeffs(4,:))
    
              ! --- Use boys function for all n = 0 : nF-1 ---
              do n = 0, nF-1
                F(n+1, p1234) = boys(n, x)
              end do
            end do
          end do
        end do
        end do
    end subroutine build_F

    function compute_boys_x(A, B, C, D, exp1, coeff1, exp2, coeff2, exp3, coeff3, exp4, coeff4) result(x)
        use ISO_C_BINDING, only: c_double
        implicit none
    
        ! Inputs
        real(c_double), dimension(3), intent(in) :: A, B, C, D
        real(c_double), dimension(:), intent(in) :: exp1, coeff1
        real(c_double), dimension(:), intent(in) :: exp2, coeff2
        real(c_double), dimension(:), intent(in) :: exp3, coeff3
        real(c_double), dimension(:), intent(in) :: exp4, coeff4
    
        ! Output
        real(c_double) :: x
    
        ! Locals
        real(c_double) :: alpha1, alpha2, alpha3, alpha4, rho1, rho2
        real(c_double) :: weight, total_weight, avg_x
        real(c_double), dimension(3) :: Pcenter, Qcenter, diff
        integer :: i1, i2, j1, j2
    
        avg_x = 0.0d0
        total_weight = 0.0d0
    
        do i1 = 1, size(exp1)
            alpha1 = exp1(i1)
            do i2 = 1, size(exp2)
                alpha2 = exp2(i2)
                rho1 = alpha1 * alpha2 / (alpha1 + alpha2)
                Pcenter = (alpha1 * A + alpha2 * B) / (alpha1 + alpha2)
    
                do j1 = 1, size(exp3)
                    alpha3 = exp3(j1)
                    do j2 = 1, size(exp4)
                        alpha4 = exp4(j2)
                        rho2 = alpha3 * alpha4 / (alpha3 + alpha4)
                        Qcenter = (alpha3 * C + alpha4 * D) / (alpha3 + alpha4)
    
                        diff = Pcenter - Qcenter
                        weight = coeff1(i1) * coeff2(i2) * coeff3(j1) * coeff4(j2) * &
                                 exp(-alpha1 - alpha2 - alpha3 - alpha4)
    
                        avg_x = avg_x + weight * sum(diff**2) / (4d0 * (rho1 + rho2))
                        total_weight = total_weight + weight
                    end do
                end do
            end do
        end do
    
        if (total_weight > 0d0) then
            x = avg_x / total_weight
        else
            x = 0.0d0
        end if
    end function compute_boys_x

end module boys

