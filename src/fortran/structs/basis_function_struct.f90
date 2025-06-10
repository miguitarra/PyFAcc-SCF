module basis_function_struct
    use iso_c_binding
    use molecule_struct
    !> Always declare everything explicitly
    implicit none

    ! Define a shell for an atom (e.g. S, P, D)
    type :: BasisShell_type
        character(len=2) :: shell_type  ! e.g., "S", "SP", "D"
        real(c_double), allocatable :: exponents(:)
        real(c_double), allocatable :: coefficients(:,:) ! (# primitives, # components in shell)
    end type
    
    ! Stores all basis shells for an atom type
    type :: AtomBasis_type
        character(len=2) :: symbol
        type(BasisShell_type), allocatable :: shells(:)
    end type
    
    ! Final structure used in SCF routines
    type :: BasisFunction_type
        integer(c_int) :: atom
        real(c_double), allocatable :: exponents(:)
        real(c_double), allocatable :: coefficients(:)
        real(c_double) :: center(3)
        integer :: ang_mom(3)
    end type
        
    type :: Shell_type
        real(c_double), allocatable :: exponents(:)
        real(c_double), allocatable :: coefficients(:)
        real(c_double) :: center(3)
        integer :: l ! e.g., S: 0, P: 1, D: 2 ...
        integer, allocatable :: bf_idx(:) ! holds the global indices of basis functions in this shell
    end type

end module basis_function_struct