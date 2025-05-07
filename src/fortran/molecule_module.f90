module molecule_module
    use iso_c_binding
    implicit none
    integer, parameter :: wp = selected_real_kind(15)

    
    type :: Molecule
        integer :: natoms
        integer :: nelectrons
        integer(c_int), allocatable :: atomic_numbers(:)
        real(c_double), allocatable :: coords(:,:) ! (3, natoms)
        integer(c_int), allocatable :: nbasis_atom(:) ! Changed to 1D array
        real(c_double), allocatable :: basis_exp(:)
        real(wp) :: Enuc ! Nuclear energy
    end type Molecule
    
contains
    
    subroutine init_molecule(mol, natoms,nelectrons, atomic_numbers, coords, nbasis_atom, basis_exp)
        type(Molecule), intent(out) :: mol
        integer(c_int), intent(in) :: natoms
        integer(c_int), intent(in) :: nelectrons
        integer(c_int), intent(in) :: atomic_numbers(:)
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: nbasis_atom(:) ! Changed to 1D array
        real(c_double), intent(in) :: basis_exp(:)
        integer :: i, j
        real(wp) :: Rij
    
        ! Initialize scalar properties
        mol%natoms = natoms
        mol%nelectrons = nelectrons
        mol%nbasis_atom = nbasis_atom
    
        ! Allocate arrays
        allocate(mol%atomic_numbers(natoms))
        allocate(mol%coords(3, natoms))
        allocate(mol%basis_exp(size(basis_exp)))
    
        ! Copy input values to molecule fields
        mol%atomic_numbers = atomic_numbers
        mol%coords = coords
        mol%basis_exp = basis_exp

        ! Nuclear energy
        mol%Enuc = 0.0_wp
        do i = 1, mol%natoms - 1
            do j = i + 1, mol%natoms
                Rij = norm2(mol%coords(:,i) - mol%coords(:,j))
                mol%Enuc = mol%Enuc + mol%atomic_numbers(i) * mol%atomic_numbers(j) / Rij
            end do
        end do
        !print *, "Nuclear energy calculated:", mol%Enuc
    end subroutine init_molecule

end module molecule_module

    