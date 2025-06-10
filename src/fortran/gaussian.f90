module gaussian
    use iso_c_binding
    use molecule_struct
    use openacc
    use constants_struct
    use utils
    implicit none
contains

    pure function primitive_norm_factor(alpha, l) result(norm)
        use iso_c_binding
        real(c_double), intent(in) :: alpha
        integer, intent(in) :: l
        real(c_double) :: norm
        real(c_double) :: pi
    
        pi = 4.0d0 * atan(1.0d0)
        norm = (2.0d0 * alpha / pi)**0.75d0 * (4.0d0 * alpha)**(0.5d0 * l) / sqrt(dble(factorial2(2*l - 1)))
    end function primitive_norm_factor

    
    pure function gaussian_overlap(alpha, beta, l) result(s)
        use iso_c_binding
        real(c_double), intent(in) :: alpha, beta
        integer, intent(in) :: l
        real(c_double) :: s
        real(c_double) :: gamma, pi
    
        pi = 4.0d0 * atan(1.0d0)
        gamma = alpha + beta
        s = (pi / gamma)**1.5d0 * factorial2(2*l - 1) / (2.0d0 * gamma)**l
    end function gaussian_overlap

    pure function gaussian_overlap3d(alpha, beta, l) result(s)
        use iso_c_binding
        real(c_double), intent(in) :: alpha, beta
        integer, intent(in) :: l
        real(c_double) :: s
        real(c_double) :: gamma, pi
    
        pi = 4.0d0 * atan(1.0d0)
        gamma = alpha + beta
        s = (pi / gamma)**1.5d0 * factorial2(2*l - 1) / (2.0d0 * gamma)**l
    end function gaussian_overlap3d

    SUBROUTINE gaussian_product(aa,bb,Ra,Rb,pp,Rp,cp)
        !$acc routine seq
        IMPLICIT NONE
        
        ! INPUT
        REAL*8, intent(in) :: aa, bb ! Gaussian exponential coefficients
        REAL*8, dimension(3), intent(in) :: Ra, Rb ! Gaussian centers

        ! INTERMEDIATE VARIABLE
        REAL*8, dimension(3) :: Rab ! Differente between Gaussian centers

        ! OUTPUT
        REAL*8, intent(out) :: pp ! Gaussian product exponential coefficient
        REAL*8, intent(out) :: cp ! Gaussian product coefficient
        REAL*8, dimension(3), intent(out) :: Rp ! Gaussian product center

        ! Compute Gaussian product exponential coefficient
        pp = aa + bb

        ! Compute difference between Gaussian centers
        Rab = Ra - Rb

        ! Compute Gaussian product coefficient
        cp = DOT_PRODUCT(Rab,Rab)
        cp =  - aa * bb / pp * cp
        cp = EXP(cp)

        ! Compute Gaussian product center
        Rp = (aa * Ra + bb * Rb) / pp

    END SUBROUTINE gaussian_product

    ! --------------------------------------------------------------------
    ! Norm of Cartesian Gaussian integrals with arbitrary angular momentum
    ! --------------------------------------------------------------------
    FUNCTION norm(ax,ay,az,aa) result(N)
        ! INPUT
        !$acc routine seq
        INTEGER, intent(in) :: ax, ay, az ! Cartesian Gaussian angular momenta projections
        REAL*8, intent(in) :: aa ! Gaussian exponential coefficient

        ! OUTPUT
        REAL*8 :: N ! Normalization factor

        integer :: l
        real*8 :: lhalf
        
        l = ax + ay + az
        lhalf = real(l, 8) / 2.0D0
        
        N = (2.0D0 * aa / PI)**(3.0D0 / 4.0D0) * (4.0D0 * aa)**lhalf
        N = N / SQRT(REAL(factorial2(2 * ax - 1) * factorial2(2 * ay - 1) * factorial2(2 * az - 1), 8))

    END FUNCTION norm


end module gaussian