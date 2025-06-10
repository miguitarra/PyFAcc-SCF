module molecule_struct
    use iso_c_binding
    use constants_struct
    implicit none

    public :: Molecule_type
    
    type :: Molecule_type
        integer :: natoms
        integer :: nelectrons
        integer :: nalpha
        integer :: nbeta
        integer(c_int), allocatable :: atomic_numbers(:)
        real(c_double), allocatable :: coords(:,:) ! (natoms, 3)
        real(wp) :: Enuc ! Nuclear energy
    end type 
end module molecule_struct