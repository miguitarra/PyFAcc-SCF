module electronic
    use constants_struct
    use nuclear
    use utils
    use openacc
    use eri_struct
    use boys
    implicit none

    real(8), parameter :: threshold = 1.0D-10
contains

    pure function ij_index(i, j) result(index)
        integer, intent(in) :: i, j
        integer :: index
        if (i >= j) then
            index = i * (i - 1) / 2 + j
        else
            index = j * (j - 1) / 2 + i
        end if
    end function ij_index


    pure function eri_indexf(mu, nu, lam, sig) result(idx)
        integer, intent(in) :: mu, nu, lam, sig
        integer :: idx
        integer :: ij, kl
    
        ij = ij_index(mu, nu)
        kl = ij_index(lam, sig)
    
        if (ij >= kl) then
            idx = ij * (ij + 1) / 2 + kl + 1  ! +1 for Fortran 1-based indexing
        else
            idx = kl * (kl + 1) / 2 + ij + 1
        end if
    end function eri_indexf

    function theta(l,l1,l2,a,b,r,g) result(t)
        !----------------------------------------
        ! Coefficient theta
        !----------------------------------------
        !
        ! Source:
        !     Handbook of Computational Chemistry
        !     David Cook
        !     Oxford University Press
        !     1998
        !
        !----------------------------------------

        ! INPUT
        integer, intent(in) :: l, l1, l2, r
        real*8, intent(in) :: a, b, g

        ! OUTPUT
        real*8 :: t

        t = 1.0D0
        t = t * f(l,l1,l2,a,b) * factorial(l) * g**(r-l)
        t = t / ( factorial(r) * factorial(l - 2*r) )

    end function theta
    
    
    function B(l,ll,r,rr,i,l1,l2,Ra,Rb,Rp,g1,l3,l4,Rc,Rd,Rq,g2,delta)
        !----------------------------------------
        ! Coefficient B
        !----------------------------------------
        !
        ! Source:
        !     Handbook of Computational Chemistry
        !     David Cook
        !     Oxford University Press
        !     1998
        !
        !----------------------------------------


        ! INPUT
        integer, intent(in) :: l, ll, r, rr, i, l1, l2, l3, l4
        real*8, intent(in) :: Ra, Rb, Rp, g1, Rc, Rd, Rq, g2, delta

        ! OUTPUT
        real*8 :: B

        B = 1.0D0
        B = B * (-1)**(l) * theta(l,l1,l2,Rp-Ra,Rp-Rb,r,g1)
        B = B * theta(ll,l3,l4,Rq-Rc,Rq-Rd,rr,g2) * (-1)**i
        B = B * (2.0D0*delta)**(2*(r+rr)) * factorial(l + ll - 2 * r - 2 * rr)
        B = B * delta**i * (Rp-Rq)**(l+ll - 2*(r+rr+i))
        B = B / ( (4.0D0*delta)**(l+ll) * factorial(i) * factorial(l+ll - 2*(r+rr+i)) )

    end function B
    
    function electronic_coeff(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,aa,bb,cc,dd,Ra,Rb,Rc,Rd) result(G)
        !$acc routine seq
        !--------------------------------------------------------------------------
        ! Compute electron-electron integral between two cartesia Gaussian function
        !--------------------------------------------------------------------------
        !
        ! Source:
        !     Handbook of Computational Chemistry
        !     David Cook
        !     Oxford University Press
        !     1998
        !
        !--------------------------------------------------------------------------
        

        ! INPUT
        integer, intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz   ! Angular momentum coefficients
        real*8, intent(in) :: aa, bb, cc, dd                                    ! Exponential Gaussian coefficients
        real*8, dimension(3), intent(in) :: Ra, Rb, Rc, Rd                      ! Gaussian centers

        ! INTERMEDIATE VARIABLES
        real*8 :: delta
        real*8 :: g1, g2                        ! Gaussian produc exponential coefficient
        real*8, dimension(3) :: Rp, Rq          ! Gaussian produc center
        real*8 :: c1, c2                        ! Gaussian product multiplicative constant
        real*8 :: BBx, BBy, BBz                 ! Temporary calls to function B
        integer :: nu                           ! Boys function index
        integer :: l, r, i, ll, rr              ! Indices loop 1
        integer :: m, s, j, mm, ss              ! Indices loop 2
        integer :: n, t, k, nn, tt              ! Indices loop 3
        real(c_double) :: vals(0:12)

        ! OUTPUT
        real*8 :: G ! Electron-electron repulsion coefficient

        G = 0.0D0

        call gaussian_product(aa,bb,Ra,Rb,g1,Rp,c1) ! Gaussian product of left gaussians
        call gaussian_product(cc,dd,Rc,Rd,g2,Rq,c2) ! Gaussian product of right gaussians

        delta = 1.0D0 / (4.0D0 * g1) + 1.0D0 / (4.0D0 * g2)

        do l = 0, ax + bx
            do r = 0, FLOOR(l / 2.0)
                do ll = 0, cx + dx
                    do rr = 0, FLOOR(ll / 2.0)
                        do i = 0, FLOOR((l+ll-2*r-2*rr) / 2.0)

                             BBx = B(l,ll,r,rr,i,ax,bx,Ra(1),Rb(1),Rp(1),g1,&
                                     cx,dx,Rc(1),Rd(1),Rq(1),g2,delta)

                             do m = 0, ay + by
                                 do s = 0, FLOOR(m / 2.0)
                                     do mm = 0, cy + dy
                                         do ss = 0, FLOOR(mm / 2.0)
                                             do j = 0, FLOOR((m+mm-2*s-2*ss) / 2.0)
                                                BBy = B(m,mm,s,ss,j,ay,by,Ra(2),Rb(2),Rp(2),&
                                                        g1,cy,dy,Rc(2),Rd(2),Rq(2),g2,delta)

                                                do n = 0, az + bz
                                                    do t = 0, FLOOR(n / 2.0)
                                                        do nn = 0, cz + dz
                                                            do tt = 0, FLOOR(nn / 2.0)
                                                                do k = 0, FLOOR((n+nn-2*t-2*tt) / 2.0)

                                                                    BBz = B(n,nn,t,tt,k,az,bz,Ra(3),Rb(3),Rp(3),g1,&
                                                                            cz,dz,Rc(3),Rd(3),Rq(3),g2,delta)

                                                                    nu = l+ll+m+mm+n+nn - 2*(r+rr+s+ss+t+tt) - (i + j + k)
                                                                    G = G + BBx * BBy * BBz * &
                                                                                            boys(nu, dot_product(Rp - Rq, Rp - Rq) * g1)

                                                                end do ! tt
                                                            end do ! nn
                                                        end do ! k
                                                    end do ! t
                                                end do ! n

                                            end do ! ss
                                        end do ! mm
                                    end do ! j
                                end do ! s
                            end do ! m

                        end do ! rr
                    end do ! ll
                end do ! i
            end do ! r
        end do ! k

        G = G * norm(ax,ay,az,aa) * norm(bx,by,bz,bb) * norm(cx,cy,cz,cc) * norm(dx,dy,dz,dd)
        G = G * c1 * c2 * 2.0D0 * PI**2 / (g1 * g2) * SQRT(PI / (g1 + g2))

    end function electronic_coeff


    subroutine compute_ERI(mol, basis_function, n_unique, ERI_flat, eri_list)
    
        ! Input
        type(Molecule_type), intent(in) :: mol
        type(BasisFunction_type), intent(in) :: basis_function(:)
        integer*8, intent(in) :: n_unique
        type(ERI_Index_type), intent(in) :: eri_list(n_unique)
    
        ! Output
        real(c_double), intent(out) :: ERI_flat(:)
    
        ! Locals
        integer :: i, p1, p2, p3, p4
        integer :: mu, nu, lam, sig
        integer :: idx_eri
        real(c_double) :: eri_val
        real(c_double) :: coeff1, coeff2, coeff3, coeff4
        real(c_double) :: exp1, exp2, exp3, exp4
    
        ! For basis function details
        type(BasisFunction_type) :: bf1, bf2, bf3, bf4
    
        do i = 1, n_unique
            mu  = eri_list(i)%mu
            nu  = eri_list(i)%nu
            lam = eri_list(i)%lam
            sig = eri_list(i)%sig
    
            bf1 = basis_function(mu)
            bf2 = basis_function(nu)
            bf3 = basis_function(lam)
            bf4 = basis_function(sig)
    
            eri_val = 0.0d0
    
            do p1 = 1, size(bf1%exponents)
                do p2 = 1, size(bf2%exponents)
                    do p3 = 1, size(bf3%exponents)
                        do p4 = 1, size(bf4%exponents)
                            exp1 = bf1%exponents(p1)
                            exp2 = bf2%exponents(p2)
                            exp3 = bf3%exponents(p3)
                            exp4 = bf4%exponents(p4)
    
                            coeff1 = bf1%coefficients(p1)
                            coeff2 = bf2%coefficients(p2)
                            coeff3 = bf3%coefficients(p3)
                            coeff4 = bf4%coefficients(p4)
    
                            eri_val = eri_val + coeff1 * coeff2 * coeff3 * coeff4 * &
                                      electronic_coeff(bf1%ang_mom(1), bf1%ang_mom(2), bf1%ang_mom(3), &
                                                       bf2%ang_mom(1), bf2%ang_mom(2), bf2%ang_mom(3), &
                                                       bf3%ang_mom(1), bf3%ang_mom(2), bf3%ang_mom(3), &
                                                       bf4%ang_mom(1), bf4%ang_mom(2), bf4%ang_mom(3), &
                                                       exp1, exp2, exp3, exp4, &
                                                       bf1%center, bf2%center, bf3%center, bf4%center)
                        end do
                    end do
                end do
            end do
    
            idx_eri = eri_indexf(mu, nu, lam, sig)
            ERI_flat(idx_eri) = eri_val
        end do
    end subroutine compute_ERI


    function compute_eri_value(mu, nu, lam, sig, mol, basis_function) result(eri_val)
        !$acc routine seq
        ! Input
        type(Molecule_type), intent(in) :: mol
        type(BasisFunction_type), intent(in) :: basis_function(:)
        
        ! Locals
        integer :: i, p1, p2, p3, p4
        integer :: mu, nu, lam, sig
        integer :: idx_eri
        real(c_double) :: eri_val
        real(c_double) :: coeff1, coeff2, coeff3, coeff4
        real(c_double) :: exp1, exp2, exp3, exp4
        real(c_double) :: product_coeffs
    
        ! For basis function details
        type(BasisFunction_type) :: bf1, bf2, bf3, bf4
    
        ! Threshold for screening
        real(c_double), parameter :: threshold = 1.0d-12
    
        bf1 = basis_function(mu)
        bf2 = basis_function(nu)
        bf3 = basis_function(lam)
        bf4 = basis_function(sig)
    
        eri_val = 0.0d0
    
        do p1 = 1, size(bf1%exponents)
            do p2 = 1, size(bf2%exponents)
                do p3 = 1, size(bf3%exponents)
                    do p4 = 1, size(bf4%exponents)
                        exp1 = bf1%exponents(p1)
                        exp2 = bf2%exponents(p2)
                        exp3 = bf3%exponents(p3)
                        exp4 = bf4%exponents(p4)
    
                        coeff1 = bf1%coefficients(p1)
                        coeff2 = bf2%coefficients(p2)
                        coeff3 = bf3%coefficients(p3)
                        coeff4 = bf4%coefficients(p4)
    
                        product_coeffs = abs(coeff1 * coeff2 * coeff3 * coeff4)
                        
                        ! === SCREENING based on product of coefficients ===
                        if (product_coeffs < threshold) cycle
    
                        eri_val = eri_val + coeff1 * coeff2 * coeff3 * coeff4 * &
                            electronic_coeff(bf1%ang_mom(1), bf1%ang_mom(2), bf1%ang_mom(3), &
                                             bf2%ang_mom(1), bf2%ang_mom(2), bf2%ang_mom(3), &
                                             bf3%ang_mom(1), bf3%ang_mom(2), bf3%ang_mom(3), &
                                             bf4%ang_mom(1), bf4%ang_mom(2), bf4%ang_mom(3), &
                                             exp1, exp2, exp3, exp4, &
                                             bf1%center, bf2%center, bf3%center, bf4%center)
                    end do
                end do
            end do
        end do
    
    end function compute_eri_value

    
end module electronic