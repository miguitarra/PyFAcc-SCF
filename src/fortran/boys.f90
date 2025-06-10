module boys
    use iso_c_binding
    use utils
    implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2021 The Numericus Group, LLC.  All rights reserved.
! You may use, distribute, and modify this code ONLY under the terms
! of the TNG IP license.
!
! You should have received a copy of the TNG IP license with this
! file.  If not, please write to: The Numericus Group, LLC,
! 2525 Arapahoe Ave. E4-431, Boulder, CO 80302.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   NAME:dboysfun12
!
!   DESC: Computes values of the Boys function for n=0,1,...  12
!         for a real argument
!
!        Input: x  --- argument, real *8, x >= 0
!
!        Output: vals  --- values of the Boys function, n = 0,1,...  12
!
!
!   DATE: MAY 2021
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine dboysfun12(x,vals)
       implicit none
       character, parameter :: name*12="dboysfun12: "
!
       real *8, parameter :: tol= 1.0d-03
       real *8, parameter :: sqrtpio2=0.886226925452758014d0
       real *8, parameter :: t(0:  11)=[ &
                        0.20000000000000000D+01,&
                        0.66666666666666663D+00,&
                        0.40000000000000002D+00,&
                        0.28571428571428570D+00,&
                        0.22222222222222221D+00,&
                        0.18181818181818182D+00,&
                        0.15384615384615385D+00,&
                        0.13333333333333333D+00,&
                        0.11764705882352941D+00,&
                        0.10526315789473684D+00,&
                        0.95238095238095233D-01,&
                        0.86956521739130432D-01]
       complex *16, parameter :: zz(1:  10)=[ &
                      (  0.64304020652330500D+01,  0.18243694739308491D+02),&
                      (  0.64304020652330500D+01, -0.18243694739308491D+02),&
                      ( -0.12572081889410178D+01,  0.14121366415342502D+02),&
                      ( -0.12572081889410178D+01, -0.14121366415342502D+02),&
                      ( -0.54103079551670268D+01,  0.10457909575828442D+02),&
                      ( -0.54103079551670268D+01, -0.10457909575828442D+02),&
                      ( -0.78720025594983341D+01,  0.69309284623985663D+01),&
                      ( -0.78720025594983341D+01, -0.69309284623985663D+01),&
                      ( -0.92069621609035313D+01,  0.34559308619699376D+01),&
                      ( -0.92069621609035313D+01, -0.34559308619699376D+01)]
       complex *16, parameter :: fact(1:  10)=[ &
                      (  0.13249210991966042D-02,  0.91787356295447745D-03),&
                      (  0.13249210991966042D-02, -0.91787356295447745D-03),&
                      (  0.55545905103006735D-01, -0.35151540664451613D+01),&
                      (  0.55545905103006735D-01,  0.35151540664451613D+01),&
                      ( -0.11456407675096416D+03,  0.19213789620924834D+03),&
                      ( -0.11456407675096416D+03, -0.19213789620924834D+03),&
                      (  0.20915556220686653D+04, -0.15825742912360638D+04),&
                      (  0.20915556220686653D+04,  0.15825742912360638D+04),&
                      ( -0.94779394228935325D+04,  0.30814443710192086D+04),&
                      ( -0.94779394228935325D+04, -0.30814443710192086D+04)]
       complex *16, parameter :: ww(1:  10)=[ &
                      ( -0.83418049867878959D-08, -0.70958810331788253D-08),&
                      ( -0.83418050437598581D-08,  0.70958810084577824D-08),&
                      (  0.82436739552884774D-07, -0.27704117936134414D-06),&
                      (  0.82436739547688584D-07,  0.27704117938414886D-06),&
                      (  0.19838416382728666D-05,  0.78321058613942770D-06),&
                      (  0.19838416382681279D-05, -0.78321058613180811D-06),&
                      ( -0.47372729839268780D-05,  0.58076919074212929D-05),&
                      ( -0.47372729839287016D-05, -0.58076919074154416D-05),&
                      ( -0.68186014282131608D-05, -0.13515261354290787D-04),&
                      ( -0.68186014282138385D-05,  0.13515261354295612D-04)]
       real *8, parameter :: rzz(1:   1)=[ &
                       -0.96321934290343840D+01]
       real *8, parameter :: rfact(1:   1)=[ &
                        0.15247844519077540D+05]
      real *8, parameter :: rww(1:   1)=[ &
                        0.18995875677635889D-04]
       real *8, intent(in)  :: x
       real *8, intent(out) :: vals(0:  12)
       real *8    :: y, yy,  rtmp
       real *8    ::   p, q, tmp
       integer *4 :: n, k
!
          y = exp(-x)
!
          if (abs(x).ge. 0.45425955121971775D+01) then
          yy      = sqrt(x)
          vals(0) = sqrtpio2*erf(yy)/yy
         yy = y/2.0d0
          do n = 1,  12
          vals(n) = ((n -0.5d0)*vals(n-1) - yy)/x 
          enddo
          return
          endif
!
          rtmp = 0
          do k = 1,   10,2
          rtmp = rtmp + ww(k)*(1.0d0 - fact(k)*y)/(x + zz(k))
          enddo
!
          tmp = 0
          do k = 1,    1
          if (abs(x + rzz(k)).ge.tol) then 
          tmp = tmp + rww(k)*(1.0d0 - rfact(k)*y)/(x + rzz(k))
          else
          q = x+rzz(k)
          p = 1.0d0 - q/2.0d0 + q**2/6.0d0 - q**3/24.0d0 + q**4/120.0d0
          tmp = tmp +rww(k)*p
          endif
          enddo
!
          vals(  12) = 2*rtmp+tmp
          yy = y/2.0d0
          do n =   11,0,-1
          vals(n) = (x*vals(n+1)+yy)*t(n)
          enddo
!
     return
     end subroutine dboysfun12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! END: dboysfun12
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\



    subroutine build_F(am1, am2, am3, am4, deriv_order, contrdepth, F, centers, alphas, coeffs)
        implicit none
        integer, intent(in) :: am1, am2, am3, am4, deriv_order, contrdepth
        ! F will be (am1+am2+am3+am4+1+deriv_order, contrdepth**4)
        real(C_DOUBLE), intent(out) :: F(:,:)
        ! Add all the parameters you need to compute x:
        ! Example: centers(4,3), alphas(4,contrdepth)
        real(8), intent(in) :: centers(4,3), alphas(4,contrdepth), coeffs(4, contrdepth)

    
        integer :: nF, Np
        integer :: n, p1, p2, p3, p4, p1234
        real(C_DOUBLE) :: x, vals(0:12)
        
        nF = am1 + am2 + am3 + am4 + 1 + deriv_order
        Np = contrdepth**4
        
        F = 0.0d0
        
        do p1 = 1, contrdepth
        do p2 = 1, contrdepth
          do p3 = 1, contrdepth
            do p4 = 1, contrdepth
              ! Flatten contraction indices for Fortran column-major ordering:
              p1234 = (p1-1)*contrdepth**3 + (p2-1)*contrdepth**2 + (p3-1)*contrdepth + p4
        
              ! --- Compute Boys function argument x for this combination ---
              x = compute_boys_x(centers(1,:), centers(2,:), centers(3,:), centers(4,:), alphas(1,:), coeffs(1,:), alphas(2,:), coeffs(2,:), alphas(3,:), coeffs(3,:), alphas(4,:), coeffs(4,:))  ! You must define this
        
              ! --- Call the Boys function routine for all n=0:12 ---
              call dboysfun12(x, vals)
        
              ! --- Fill the F array for this contraction combination ---
              do n = 0, nF-1
                F(n+1, p1234) = vals(n)
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

