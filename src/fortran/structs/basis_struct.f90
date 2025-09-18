module basis_struct
    use iso_c_binding
    use constants_struct
    !> Always declare everything explicitly
    implicit none

    ! Define shells
    type :: Shell_type
        integer :: atom_index ! connect to the atom
        integer :: l_num ! e.g., "S" = 0, "P"= 1, "D" = 2
        integer :: num_prim ! for contraction purposes
        real(wp), allocatable :: exponents(:)
        real(wp), allocatable :: coefficients (:)
        real(wp) :: x, y, z 
    end type

    type :: Atom_type
        integer :: z_num ! atomic number
        real(wp) :: x, y, z ! center coordinates
        integer :: n_shells ! number of shells, probably not necessary but just in case :)
    end type

    type :: Molecule_type
        integer :: n_atoms ! number of atoms in the system
        type(Atom_type), allocatable :: atoms(:) ! array of atoms
        integer :: charge ! charge of the system
        integer :: nalpha ! number of alpha electrons
        integer :: nbeta ! number of beta electrons
        integer :: nelectrons ! total number of electrons
        !integer :: multiplicity ! multiplicity of the system
        real(wp) :: Enuc ! nuclear energy, calculated in the molecule.f90 module
    end type

end module basis_struct