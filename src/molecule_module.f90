module molecule_module
    use iso_c_binding
    implicit none
    
    type :: Molecule
        integer :: natoms
        integer(c_int), allocatable :: atomic_numbers(:)
        real(c_double), allocatable :: coords(:,:) ! (3, natoms)
        integer(c_int), allocatable :: nbasis_atom(:) ! Changed to 1D array
        real(c_double), allocatable :: basis_exp(:)
    end type Molecule
    
contains
    
    subroutine init_molecule(mol, natoms, atomic_numbers, coords, nbasis_atom, basis_exp)
        type(Molecule), intent(out) :: mol
        integer(c_int), intent(in) :: natoms
        integer(c_int), intent(in) :: atomic_numbers(:)
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: nbasis_atom(:) ! Changed to 1D array
        real(c_double), intent(in) :: basis_exp(:)
    
        ! Initialize scalar properties
        mol%natoms = natoms
        mol%nbasis_atom = nbasis_atom
    
        ! Allocate arrays
        allocate(mol%atomic_numbers(natoms))
        allocate(mol%coords(3, natoms))
        allocate(mol%basis_exp(size(basis_exp)))
    
        ! Copy input values to molecule fields
        mol%atomic_numbers = atomic_numbers
        mol%coords = coords
        mol%basis_exp = basis_exp
    end subroutine init_molecule

end module molecule_module

    