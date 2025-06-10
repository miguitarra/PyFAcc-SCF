module basis
    use iso_c_binding
    use molecule_struct
    use basis_function_struct
    use utils
    use gaussian
    implicit none

contains

    subroutine build_basis_functions(atom_symbols, coords, natoms, basis_file, basis_functions, shell)
        character(len=*), dimension(natoms), intent(in) :: atom_symbols
        real(c_double), dimension(3, natoms), intent(in) :: coords
        integer(c_int), intent(in) :: natoms
        character(len=*), intent(in) :: basis_file
        type(BasisFunction_type), allocatable, intent(out) :: basis_functions(:)
        type(Shell_type), allocatable, intent(out) :: shell(:)
        
        type(AtomBasis_type), allocatable :: atom_bases(:)
        integer :: i, j, k, nbf, count, p, d, count_shell, nshells
        character(len=2) :: atom
        integer :: ncomponents

        !print *, '>> Entering build_basis_functions'
        call parse_basis_file(trim(basis_file), atom_bases)
        !print *, '>> Parsed basis file'

        ! Count total basis functions
        nbf = 0
        nshells = 0
        do i = 1, natoms
            atom = trim(atom_symbols(i))
            !print *, '>> Processing atom: ', atom
            do j = 1, size(atom_bases)
                if (trim(atom_bases(j)%symbol) == atom) then
                    !print *, '   Found matching atom in basis file: ', atom
                    do k = 1, size(atom_bases(j)%shells)
                        select case (trim(atom_bases(j)%shells(k)%shell_type))
                        case ('SP')
                            nbf = nbf + 4  ! 1 S + 3 P functions
                            nshells = nshells + 2
                        case ('D')
                            nbf = nbf + 6  ! 5 D functions
                            nshells = nshells +1
                        case default
                            nbf = nbf + 1
                            nshells = nshells +1
                        end select
                    end do
                    exit
                end if
            end do
        end do

        !print *, '>> Total basis functions to allocate: ', nbf
        allocate(basis_functions(nbf))
        allocate(shell(nshells))
        count = 0
        count_shell = 0

        do i = 1, natoms
            atom = trim(atom_symbols(i))
            do j = 1, size(atom_bases)
                if (trim(atom_bases(j)%symbol) == atom) then
                    do k = 1, size(atom_bases(j)%shells)
                        ncomponents = size(atom_bases(j)%shells(k)%coefficients, 2)
                        !print *, '   Shell type: ', trim(atom_bases(j)%shells(k)%shell_type), ' with ', ncomponents, ' components'

                        select case (trim(atom_bases(j)%shells(k)%shell_type))
                        case ('SP')
                            ! ----- S section -----
                            ! ------ Basis function -----
                            count = count + 1
                            count_shell = count_shell + 1
                            !print *, '   -> Allocating S from SP for basis function ', count
                            basis_functions(count)%atom = i
                            basis_functions(count)%center = coords(:,i)
                            allocate(basis_functions(count)%exponents(size(atom_bases(j)%shells(k)%exponents)))
                            allocate(basis_functions(count)%coefficients(size(atom_bases(j)%shells(k)%exponents)))
                            basis_functions(count)%exponents = atom_bases(j)%shells(k)%exponents
                            basis_functions(count)%coefficients = atom_bases(j)%shells(k)%coefficients(:, 1)
                            basis_functions(count)%ang_mom = [0,0,0]
                            
                            !------ Shell -----
                            shell(count_shell)%center = coords(:,i)
                            allocate(shell(count_shell)%exponents(size(atom_bases(j)%shells(k)%exponents)))
                            allocate(shell(count_shell)%coefficients(size(atom_bases(j)%shells(k)%exponents)))
                            shell(count_shell)%exponents = atom_bases(j)%shells(k)%exponents
                            shell(count_shell)%coefficients = atom_bases(j)%shells(k)%coefficients(:, 1)
                            shell(count_shell)%l = 0
                            allocate(shell(count_shell)%bf_idx(1))
                            shell(count_shell)%bf_idx(1) = count
                            call normalize_shell(shell(count_shell)%l, shell(count_shell)%exponents, shell(count_shell)%coefficients)

                            ! ----- P section -----
                            ! ------ Shell -----
                            count_shell = count_shell + 1
                            shell(count_shell)%center = coords(:,i)
                            allocate(shell(count_shell)%exponents(size(atom_bases(j)%shells(k)%exponents)))
                            allocate(shell(count_shell)%coefficients(size(atom_bases(j)%shells(k)%exponents)))
                            shell(count_shell)%exponents = atom_bases(j)%shells(k)%exponents
                            shell(count_shell)%coefficients = atom_bases(j)%shells(k)%coefficients(:, 1)
                            shell(count_shell)%l = 1
                            call normalize_shell(shell(count_shell)%l, shell(count_shell)%exponents, shell(count_shell)%coefficients)
                            
                            
                            allocate(shell(count_shell)%bf_idx(3))

                            ! ------ Basis function -----
                            do p = 1, 3
                                count = count + 1
                                basis_functions(count)%atom = i
                                basis_functions(count)%center = coords(:,i)
                                allocate(basis_functions(count)%exponents(size(atom_bases(j)%shells(k)%exponents)))
                                allocate(basis_functions(count)%coefficients(size(atom_bases(j)%shells(k)%exponents)))
                                basis_functions(count)%exponents = atom_bases(j)%shells(k)%exponents
                                basis_functions(count)%coefficients = atom_bases(j)%shells(k)%coefficients(:, 2)
                                
                                basis_functions(count)%ang_mom = 0
                                basis_functions(count)%ang_mom(p) = 1

                                shell(count_shell)%bf_idx(p) = count
                                
                            end do
                            
                        ! ----- D section -----
                        case ('D')
                            ! ------ Shell -----
                            count_shell = count_shell + 1
                            shell(count_shell)%center = coords(:,i)
                            allocate(shell(count_shell)%exponents(size(atom_bases(j)%shells(k)%exponents)))
                            allocate(shell(count_shell)%coefficients(size(atom_bases(j)%shells(k)%exponents)))
                            shell(count_shell)%exponents = atom_bases(j)%shells(k)%exponents
                            shell(count_shell)%coefficients = atom_bases(j)%shells(k)%coefficients(:, 1)
                            shell(count_shell)%l = 2
                            call normalize_shell(shell(count_shell)%l, shell(count_shell)%exponents, shell(count_shell)%coefficients)
                            allocate(shell(count_shell)%bf_idx(6))
                            
                            do d = 1, 6
                                count = count + 1
    
                                shell(count_shell)%bf_idx(d) = count
                                
                                ! ------ Basis function -----
                                basis_functions(count)%atom = i
                                basis_functions(count)%center = coords(:,i)
                                allocate(basis_functions(count)%exponents(size(atom_bases(j)%shells(k)%exponents)))
                                allocate(basis_functions(count)%coefficients(size(atom_bases(j)%shells(k)%exponents)))
                                basis_functions(count)%exponents = atom_bases(j)%shells(k)%exponents
                                basis_functions(count)%coefficients = atom_bases(j)%shells(k)%coefficients(:, 2)
                                select case(d)
                                case(1); basis_functions(count)%ang_mom = [2, 0, 0]  ! xx
                                case(2); basis_functions(count)%ang_mom = [0, 2, 0]  ! yy
                                case(3); basis_functions(count)%ang_mom = [0, 0, 2]  ! zz
                                case(4); basis_functions(count)%ang_mom = [1, 1, 0]  ! xy
                                case(5); basis_functions(count)%ang_mom = [1, 0, 1]  ! xz
                                case(6); basis_functions(count)%ang_mom = [0, 1, 1]  ! yz
                                end select
                                
                            end do

                        case default ! -> simple S shell
                            count = count + 1
                            count_shell = count_shell + 1

                            ! ------ Basis function -----
                            basis_functions(count)%atom = i
                            basis_functions(count)%center = coords(:,i)
                            allocate(basis_functions(count)%exponents(size(atom_bases(j)%shells(k)%exponents)))
                            allocate(basis_functions(count)%coefficients(size(atom_bases(j)%shells(k)%exponents)))
                            basis_functions(count)%exponents = atom_bases(j)%shells(k)%exponents
                            basis_functions(count)%coefficients = atom_bases(j)%shells(k)%coefficients(:, 1)
                            basis_functions(count)%ang_mom = [0,0,0]

                            ! ------ Shell -----
                            shell(count_shell)%center = coords(:,i)
                            allocate(shell(count_shell)%exponents(size(atom_bases(j)%shells(k)%exponents)))
                            allocate(shell(count_shell)%coefficients(size(atom_bases(j)%shells(k)%exponents)))
                            shell(count_shell)%exponents = atom_bases(j)%shells(k)%exponents
                            shell(count_shell)%coefficients = atom_bases(j)%shells(k)%coefficients(:, 1)
                            shell(count_shell)%l = 0
                            call normalize_shell(shell(count_shell)%l, shell(count_shell)%exponents, shell(count_shell)%coefficients)
                            
                            allocate(shell(count_shell)%bf_idx(1))
                            shell(count_shell)%bf_idx = count 
                            
                            
                        end select
                    end do
                    exit
                end if
            end do
        end do
        !print *, '>> Finished building basis functions'
    end subroutine build_basis_functions

    subroutine parse_basis_file(filename, atom_bases)
        character(len=*), intent(in) :: filename
        type(AtomBasis_type), allocatable, intent(out) :: atom_bases(:)

        character(len=256) :: line
        character(len=2) :: current_atom, shell_type
        integer :: iunit, ierr, nprim, iprim
        type(AtomBasis_type), allocatable :: temp(:)
        type(BasisShell_type), allocatable :: shells(:)
        real(c_double), allocatable :: exps(:), coeffs(:,:)
        logical :: reading_shell
        integer :: shell_count, atom_count
        integer :: dummy
        real(c_double) :: scale

        !print *, '>> Opening basis file: ', trim(filename)
        open(newunit=iunit, file=filename, status='old', action='read', iostat=ierr)
        if (ierr /= 0) stop 'Error opening basis set file'

        allocate(temp(100))
        atom_count = 0
        reading_shell = .false.

        do
            read(iunit, '(A)', iostat=ierr) line
            if (ierr /= 0) exit
            line = adjustl(line)
            if (line == '' .or. line(1:1) == '!') cycle

            ! Atom header
            if (scan(line, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ') == 1) then
                read(line, *, iostat=ierr) current_atom, dummy
                if (ierr == 0 .and. dummy == 0) then
                    if (reading_shell) then
                        temp(atom_count)%shells = shells(:shell_count)
                        deallocate(shells)
                        reading_shell = .false.
                    end if
                    atom_count = atom_count + 1
                    !print *, '>> Found atom entry: ', trim(current_atom)
                    temp(atom_count)%symbol = trim(current_atom)
                    allocate(shells(20))
                    shell_count = 0
                    reading_shell = .true.
                    cycle
            end if
            end if

            ! Shell block
            if (index(line, 'SP') == 1) then
                read(line, *, iostat=ierr) shell_type, nprim, scale
                !print *, '   -> SP shell: ', nprim, ' primitives'
                shell_count = shell_count + 1
                allocate(exps(nprim))
                allocate(coeffs(nprim, 2))
                do iprim = 1, nprim
                    read(iunit, *, iostat=ierr) exps(iprim), coeffs(iprim, 1), coeffs(iprim, 2)
                    if (ierr /= 0) stop 'Error reading SP primitive'
                end do
            else if (index(line, 'S') == 1 .or. index(line, 'D') == 1) then
                read(line, *, iostat=ierr) shell_type, nprim, scale
                !print *, '   -> ', trim(shell_type), ' shell: ', nprim, ' primitives'
                shell_count = shell_count + 1
                allocate(exps(nprim))
                allocate(coeffs(nprim, 1))
                do iprim = 1, nprim
                    read(iunit, *, iostat=ierr) exps(iprim), coeffs(iprim, 1)
                    if (ierr /= 0) stop 'Error reading primitive'
                end do
            else if (trim(line) == '****') then
                if (reading_shell) then
                    !print *, '>> Finished atom block: ', trim(temp(atom_count)%symbol)
                    temp(atom_count)%shells = shells(:shell_count)
                    deallocate(shells)
                    reading_shell = .false.
                end if
                cycle
            else
                cycle
            end if

            shells(shell_count)%shell_type = trim(shell_type)
            shells(shell_count)%exponents = exps
            shells(shell_count)%coefficients = coeffs
            deallocate(exps, coeffs)
        end do

        close(iunit)
        if (reading_shell) then
            !print *, '>> Final atom block completed: ', trim(temp(atom_count)%symbol)
            temp(atom_count)%shells = shells(:shell_count)
            deallocate(shells)
        end if

        atom_bases = temp(:atom_count)
        deallocate(temp)
        !print *, '>> Done parsing basis file'
    end subroutine parse_basis_file

    subroutine normalize_basis_function(ang_mom, exponents, coefficients)
        integer, intent(in) :: ang_mom
        real(c_double), dimension(:), intent(in) :: exponents
        real(c_double), dimension(:), intent(inout) :: coefficients
    
        integer :: i, j
        real(c_double) :: norm_const, overlap
        integer :: l
        real(c_double) :: alpha_i, alpha_j, c_i, c_j
    
        l = ang_mom
        ! Normalize each primitive
        do i = 1, size(exponents)
            coefficients(i) = coefficients(i) * primitive_norm_factor(exponents(i), l)
        end do
    
        ! Normalize the full contracted basis function
        overlap = 0.0
        do i = 1, size(exponents)
            do j = 1, size(exponents)
                alpha_i = exponents(i)
                alpha_j = exponents(j)
                c_i = coefficients(i)
                c_j = coefficients(j)
    
                overlap = overlap + c_i * c_j * gaussian_overlap3d(alpha_i, alpha_j, l)
            end do
        end do
    
        norm_const = 1.0d0 / sqrt(overlap)
        coefficients = coefficients * norm_const
    
    end subroutine normalize_basis_function

    subroutine normalize_shell(l, exponents, coefficients)
        ! Normalizes a contracted Gaussian shell (all primitives and the contraction)
        integer, intent(in) :: l
        real(c_double), intent(in) :: exponents(:)
        real(c_double), intent(inout) :: coefficients(:)
    
        integer :: i, j
        real(c_double) :: norm, norm2, prefac, alpha_i, alpha_j
        real(c_double), parameter :: pi = 3.14159265358979323846_c_double
    
        ! Primitive normalization factor for general angular momentum l
        do i = 1, size(exponents)
            alpha_i = exponents(i)
            prefac = (2.0_c_double**(2*l + 1) * alpha_i**(l + 1.5_c_double)) / (pi**1.5_c_double * double_factorial(2*l - 1))
            norm = sqrt(prefac)
            coefficients(i) = coefficients(i) * norm
        end do
    
        ! Contracted normalization (normalize the linear combination)
        norm2 = 0.0_c_double
        do i = 1, size(exponents)
            alpha_i = exponents(i)
            do j = 1, size(exponents)
                alpha_j = exponents(j)
                norm2 = norm2 + coefficients(i) * coefficients(j) * overlap_ss(alpha_i, alpha_j, l)
            end do
        end do
        norm2 = 1.0_c_double / sqrt(norm2)
        coefficients = coefficients * norm2
    
    contains
    
        ! Computes the overlap integral between two normalized primitives with angular momentum l
        real(c_double) function overlap_ss(alpha, beta, l)
            real(c_double), intent(in) :: alpha, beta
            integer, intent(in) :: l
            real(c_double), parameter :: pi = 3.14159265358979323846_c_double
            overlap_ss = (pi / (alpha + beta))**1.5_c_double * ( (2*sqrt(alpha*beta)/(alpha+beta))**l )
        end function overlap_ss
    
        ! Computes the double factorial: (2l-1)!! for integer l >= 0
        integer function double_factorial(n)
            integer, intent(in) :: n
            integer :: k
            if (n <= 0) then
                double_factorial = 1
            else
                double_factorial = 1
                do k = n, 1, -2
                    double_factorial = double_factorial * k
                end do
            end if
        end function double_factorial
    
    end subroutine normalize_shell


end module basis
