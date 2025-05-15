module scf_interface
    use iso_c_binding
    use molecule_module
    use slater
    use constants_struct
    use molecule_struct
    use scf_module
    use openacc
    implicit none

contains

    subroutine run_scf_interface_py(natoms, nelectrons, nalpha, nbeta,  atomic_numbers, coords, eps_scf, max_iter, final_energy, nbasis_atom, basis_exp, basis_set) bind(C, name="run_scf_interface_py")

        implicit none
        integer(c_int), intent(in) :: natoms
        integer(c_int), intent(in) :: nelectrons
        integer(c_int), intent(in) :: nalpha
        integer(c_int), intent(in) :: nbeta
        integer(c_int), intent(in) :: atomic_numbers(natoms)
        real(c_double), intent(in) :: coords(3,natoms)
        real(c_double), intent(in) :: eps_scf
        integer(c_int), intent(in) :: max_iter
        real(c_double), intent(out) :: final_energy
        integer(c_int), intent(in) :: nbasis_atom(natoms)
        real(c_double), intent(in) :: basis_exp(sum(nbasis_atom))
        integer(c_int), intent(in) :: basis_set

        type(Molecule) :: mol
        type(BasisFunction), allocatable :: basis_functions(:)
        integer :: i, nbf, j, atom_idx
        real(wp), allocatable :: alpha(:), coeff(:)

        ! --- Check where it is being run ---
        integer :: dev_type
        dev_type = acc_get_device_type()

        select case (dev_type)
        case (acc_device_nvidia)
            !print *, "OpenACC: Running on NVIDIA GPU."
        case (acc_device_host)
            !Only let know if it is running on host
            print *, "OpenACC: Running on host (CPU)."
        case default
            !print *, "OpenACC: Unknown or unsupported device type."
        end select
        ! ---------------------------

        ! Initialize molecule
        call init_molecule(mol, natoms, nelectrons, nalpha, nbeta, atomic_numbers, coords, nbasis_atom, basis_exp)
        !print *, "Molecule initialized successfully."

        ! Allocate room for all basis functions
        allocate(basis_functions(sum(nbasis_atom)))
        !print *, "Basis functions allocated successfully."
        nbf = 0
        atom_idx = 1  ! Start index for basis_exp

        ! Fill basis functions
        ! First do per atom
        ! Allocate temporary arrays
        allocate(alpha(basis_set), coeff(basis_set))
        do i = 1, natoms
            do j = 1, nbasis_atom(i)
                nbf = nbf + 1

                ! Expand Slater-type orbital into gaussians Gaussians
                call expand_slater(basis_set, basis_exp(atom_idx), alpha, coeff)
                !print *, "expand_slater called for atom", i, "basis", j
                !print *, "alpha = ", alpha
                !print *, "coeff = ", coeff
                !print *, "sum(coeff**2) = ", sum(coeff**2)
        
                ! Allocate and copy into basis_functions(nbf)
                allocate(basis_functions(nbf)%exponents(basis_set))
                allocate(basis_functions(nbf)%coefficients(basis_set))
                basis_functions(nbf)%exponents = alpha
                basis_functions(nbf)%coefficients = coeff
                basis_functions(nbf)%center = mol%coords(:, i)
                basis_functions(nbf)%atom = i

        
                atom_idx = atom_idx + 1
            end do
        end do
        deallocate(alpha, coeff)
        
        !print *, "Basis functions filled successfully."

        ! Run SCF
        call run_scf(mol, basis_functions, eps_scf, max_iter, final_energy)
        !print *, "SCF calculation completed successfully."

    end subroutine run_scf_interface_py

end module scf_interface