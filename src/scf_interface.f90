subroutine run_scf_interface_py(natoms, atomic_numbers, coords, eps_scf, max_iter, final_energy, nbasis_atom, basis_exp) bind(C, name="run_scf_interface_py")
    use scf_module
    use molecule_module
    use slater
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: natoms
    integer(c_int), intent(in) :: atomic_numbers(natoms)
    real(c_double), intent(in) :: coords(3,natoms)
    real(c_double), intent(in) :: eps_scf
    integer(c_int), intent(in) :: max_iter
    real(c_double), intent(out) :: final_energy
    integer(c_int), intent(in) :: nbasis_atom(natoms)
    real(c_double), intent(in) :: basis_exp(sum(nbasis_atom))


    type(Molecule) :: mol
    type(BasisFunction), allocatable :: basis_functions(:)
    integer :: i, nbf, j, atom_idx, ng
    real(wp), allocatable :: alpha(:), coeff(:)

    ! Initialize molecule
    call init_molecule(mol, natoms, atomic_numbers, coords, nbasis_atom, basis_exp)

    ! Allocate room for all basis functions
    allocate(basis_functions(natoms * sum(nbasis_atom)))
    nbf = 0
    atom_idx = 1  ! Start index for basis_exp

    ! Fill basis functions
    ! First do per atom
    do i = 1, natoms
        ng = nbasis_atom(i)
    
        do j = 1, ng
            ! Now per basis function
            nbf = nbf + 1
    
            ! Allocate temporary arrays for alpha and coeff
            allocate(alpha(1), coeff(1))
    
            ! Call expand_slater for each basis function
            call expand_slater(1, basis_exp(atom_idx), alpha, coeff)
    
            ! Assign values to the basis_functions array
            basis_functions(nbf)%atom = mol%atomic_numbers(i)
            basis_functions(nbf)%exponents = alpha(1)
            basis_functions(nbf)%coefficients = coeff(1)
    
            ! Deallocate temporary arrays
            deallocate(alpha, coeff)
    
            atom_idx = atom_idx + 1  ! Move to the next basis function
        end do
    end do

    ! Run SCF
    call run_scf(mol, basis_functions, eps_scf, max_iter, final_energy)

end subroutine run_scf_interface_py
