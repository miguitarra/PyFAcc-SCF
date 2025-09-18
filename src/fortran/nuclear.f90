module nuclear
    use constants_struct
    use gaussian
    use overlap
    use basis_struct
    use openacc
    use cartesian
    use la
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

    subroutine compute_nuclear_shell_pair(xA, xB, yA, yB, zA, zB, lA, lB, startA, endA, startB, endB, exponents, coefficients, atom_x, atom_y, atom_z, Zn, V_block)
        !$acc routine seq
        ! Inputs
        integer, intent(in) :: lA, lB, startA, endA, startB, endB, Zn
        real(wp), intent(in) :: xA, xB, yA, yB, zA, zB, atom_x, atom_y, atom_z
        real(wp), intent(in) :: exponents(:), coefficients(:)
    
        ! Output
        real(wp), intent(out) :: V_block(:,:)   ! (nA, nB)
    
        ! Local variables
        integer :: nA, nB
        integer :: muA, muB
        integer :: axA, ayA, azA, axB, ayB, azB
        integer :: kA, kB
        real(wp) :: tmp
        real(wp), dimension(3) :: RA, RB, Rn
    
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
        Rn = [atom_x, atom_y, atom_z]
    
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
                              * nuclear_coeff(axA, ayA, azA, axB, ayB, azB, &
                                                 exponents(kA), exponents(kB), RA, RB, Rn, Zn)
                    end do
                end do
                V_block(muA, muB) = tmp
            end do
        end do
    
    end subroutine compute_nuclear_shell_pair

! ------------------------------------------------------------------------


    subroutine V_nuclear(mol, shells, V)
        ! input
        type(Molecule_type), intent(in) :: mol
        type(Shell_type), intent(in) :: shells(:)
        real(wp), intent(out) :: V(:,:)
        
        ! local
        integer :: n_shells, sA, sB, nA, nB, muA, muB, gA, gB, idx, i, j, n, Nn
        integer :: prim_index, coeffA_start, coeffA_end, coeffB_start, coeffB_end
        integer, allocatable :: shell_start(:)
        integer, allocatable :: prim_start(:)
        integer :: total_prim
        real(wp), allocatable :: V_block(:,:)
        real(wp), allocatable :: shell_x(:), shell_y(:), shell_z(:), shell_coeff(:), shell_exp(:), atom_x(:), atom_y(:), atom_z(:)
        integer, allocatable :: Z_atom(:)
        integer, allocatable :: shell_prim(:), shell_l(:)
        integer :: n_upper


        real(wp), allocatable :: Vtmp(:,:,:)
        

        
        n_shells = size(shells)
        Nn = mol%n_atoms
        allocate(shell_start(n_shells+1))
        allocate(shell_prim(n_shells), shell_l(n_shells), shell_x(n_shells), shell_y(n_shells), shell_z(n_shells))
    
        ! compute total number of primitives and prim_start array
        total_prim = 0
        do sA = 1, n_shells
            total_prim = total_prim + shells(sA)%num_prim
        end do
        allocate(shell_coeff(total_prim), shell_exp(total_prim))
        allocate(prim_start(n_shells+1))

    
        ! fill arrays to make it OpenACC compatible - shell side
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

        ! fill arrays to make it OpenACC compatible - mol/atoms
        print *, "got here!"
        atom_x = mol%atoms(:)%x
        print *, atom_x
        atom_y = mol%atoms(:)%y
        atom_z = mol%atoms(:)%z
        Z_atom = mol%atoms(:)%z_num


        V(:,:) = 0.0_wp
        allocate(Vtmp(Nn, size(V, 1), size(V, 2)))
        Vtmp = 0.0_wp
        n_upper = n_shells * (n_shells + 1) / 2
    
        
        !$acc data copyin(shell_x, shell_y, shell_z, shell_l, shell_coeff, shell_prim, shell_exp, shell_start, prim_start, atom_x, atom_y, atom_z, Z_atom) create(Vtmp) async(3)
        !$acc parallel loop collapse(2) private(sA, sB, nA, nB, coeffA_start, coeffA_end, coeffB_start, coeffB_end, muA, muB, gA, gB) default(present) async(3)
        do n = 1, Nn
            do idx = 1, n_upper
                call idx_upper_triangular(n_shells, idx, sA, sB)
                nA = (shell_l(sA) + 1)*(shell_l(sA) + 2)/2
                nB = (shell_l(sB) + 1)*(shell_l(sB) + 2)/2
        
                coeffA_start = prim_start(sA)
                coeffA_end   = prim_start(sA+1) - 1
                coeffB_start = prim_start(sB)
                coeffB_end   = prim_start(sB+1) - 1
        
                allocate(V_block(nA, nB))
        
                call compute_nuclear_shell_pair(shell_x(sA), shell_x(sB), &
                                               shell_y(sA), shell_y(sB), &
                                               shell_z(sA), shell_z(sB), &
                                               shell_l(sA), shell_l(sB), &
                                               coeffA_start, coeffA_end, &
                                               coeffB_start, coeffB_end, &
                                               shell_exp, shell_coeff, atom_x(n), atom_y(n), atom_z(n), Z_atom(n), V_block)
        
                do muA = 1, nA
                    gA = shell_start(sA) + muA - 1
                    do muB = 1, nB
                        gB = shell_start(sB) + muB - 1
        
                        Vtmp(n,gA, gB) = V_block(muA, muB)
                    end do
                end do
        
                deallocate(V_block)
        
            end do
        end do
        !$acc end parallel loop
        !$acc update host(Vtmp)

        V = 0.0_wp
        do n = 1, Nn
            V = V + Vtmp(n, :, :)
        end do
    
        do i = 1, size(V, 1)
            do j = i+1, size(V, 1)
                V(j, i) = V(i, j)
            end do
        end do

        !$acc end data
        
    end subroutine V_nuclear

    
end module nuclear
