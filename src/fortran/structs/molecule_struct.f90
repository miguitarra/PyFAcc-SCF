module molecule_struct
    use iso_c_binding
    use constants_struct
    implicit none

    public :: Molecule
    
    type :: Molecule
        integer :: natoms
        integer :: nelectrons
        integer :: nalpha
        integer :: nbeta
        integer(c_int), allocatable :: atomic_numbers(:)
        real(c_double), allocatable :: coords(:,:) ! (3, natoms)
        integer(c_int), allocatable :: nbasis_atom(:) ! Changed to 1D array
        real(c_double), allocatable :: basis_exp(:)
        real(wp) :: Enuc ! Nuclear energy
    end type Molecule
end module molecule_struct