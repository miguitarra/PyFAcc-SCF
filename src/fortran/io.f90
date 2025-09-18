module io
    use iso_c_binding
    use constants_struct
    use iso_fortran_env
    use basis_struct
    implicit none
contains
    subroutine print_matrix(matrix, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: matrix(n, n)
        integer :: i, j
        
        do i = 1, n
        write(*,'(100(f12.6,1x))') (matrix(i, j), j = 1, n)
        end do
    end subroutine print_matrix

    subroutine print_vector(vec, n)
        implicit none
        integer, intent(in) :: n
        real(c_double), intent(in) :: vec(n)
        integer :: i
    
        write(*,'(100(f12.6,1x))') (vec(i), i = 1, n)
    end subroutine print_vector

    subroutine print_molecule(mol, shells)
        use iso_fortran_env, only: wp => real64
        implicit none
    
        ! Arguments
        type(Molecule_type), intent(in) :: mol
        type(Shell_type), allocatable, intent(in) :: shells(:)
    
        ! Locals
        integer :: i, j, atom_idx
        type(Atom_type) :: atom_ref
        character(len=1) :: l_symbol
    
        print *, "=============================="
        print *, " Molecule Summary"
        print *, "=============================="
        print "(A,I4)", " Number of atoms: ", mol%n_atoms
        print "(A,I4)", " Charge: ", mol%charge
        print "(A,I4)", " Electrons (total): ", mol%nelectrons
        print "(A,I4)", " Alpha electrons: ", mol%nalpha
        print "(A,I4)", " Beta electrons : ", mol%nbeta
        print "(A,F12.6)", " Nuclear repulsion energy: ", mol%Enuc
        print *, "------------------------------"
    
        ! Print atoms
        do i = 1, mol%n_atoms
            print "(A,I4,A,I4)", " Atom ", i, " (Z = ", mol%atoms(i)%z_num, ")"
            print "(A,3F12.6)", "   Coordinates (Bohr): ", &
                                mol%atoms(i)%x, mol%atoms(i)%y, mol%atoms(i)%z
            print "(A,I4)", "   Number of shells: ", mol%atoms(i)%n_shells
        end do
    
        print *, "=============================="
        print *, " Shells"
        print *, "=============================="
    
        if (allocated(shells)) then
            do i = 1, size(shells)
                atom_idx = shells(i)%atom_index
                atom_ref = mol%atoms(atom_idx)
    
                select case (shells(i)%l_num)
                case (0); l_symbol = "s"
                case (1); l_symbol = "p"
                case (2); l_symbol = "d"
                case default; l_symbol = "?"
                end select
    
                print "(A,I4)", " Shell ", i
                print "(A,I4)", "   Atom index: ", atom_idx
                print "(A,I4)", "   Z number: ", atom_ref%z_num
                print "(A,A)",  "   Angular momentum: ", l_symbol
                print "(A,I4)", "   Number of primitives: ", shells(i)%num_prim

                print "(A,3F12.6)", "   Coordinates (Bohr): ", &
                                shells(i)%x, shells(i)%y, shells(i)%z
    
                if (allocated(shells(i)%exponents)) then
                    print *, "   Exponents:"
                    do j = 1, size(shells(i)%exponents)
                        print "(A,I4,F12.6)", "     ", j, shells(i)%exponents(j)
                    end do
                end if
    
                if (allocated(shells(i)%coefficients)) then
                    print *, "   Coefficients:"
                    do j = 1, size(shells(i)%coefficients)
                        print "(A,I4,F12.6)", "     ", j, shells(i)%coefficients(j)
                    end do
                end if
    
                print *, "------------------------------"
            end do
        else
            print *, "   (No shells allocated)"
        end if
    
        print *, "=============================="
        print *, " End of Molecule"
        print *, "=============================="
    end subroutine print_molecule



    function itoa(i) result(str)
        integer, intent(in) :: i
        character(len=12) :: str
        
        write(str, '(I0)') i
    end function itoa

end module io