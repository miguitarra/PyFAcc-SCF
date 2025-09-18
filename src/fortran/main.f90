module main
    use iso_c_binding
    use io
    use basis_struct
    use constants_struct
    use basis
    use molecule
    use openacc
    use rscf
    !use uscf
    implicit none

contains

    subroutine run_scf_interface_py(natoms, nelectrons, charge, nalpha, nbeta, species, coords, eps_scf, max_iter, final_energy, basis_set, len_basis_set, z_num) bind(C, name="run_scf_interface_py")

        integer(c_int32_t), intent(in) :: natoms, nelectrons, charge, nalpha, nbeta, max_iter
        character(*), intent(in) :: basis_set
        integer(c_int), intent(in) :: len_basis_set
        character(len=2), dimension(natoms), intent(in) :: species
        real(c_double), intent(in) :: coords(natoms, 3)
        integer(c_int32_t), intent(in) :: z_num(natoms)   
        real(c_double), intent(in) :: eps_scf
        real(c_double), intent(out) :: final_energy
        
        ! local variables
        type(Molecule_type) :: mol
        type(Shell_type), allocatable :: shells(:)
        character(len=200) :: basis_file
        integer :: dev_type
    
        ! define the file path
        basis_file = './data/basis_set/' // trim(basis_set(1:len_basis_set)) // '.txt'
        
        ! check if running on GPU
        dev_type = acc_get_device_type()
        select case (dev_type)
        case (acc_device_nvidia)
            !print *, "OpenACC: Running on NVIDIA GPU."
        case (acc_device_host)
            print *, "OpenACC: Running on host (CPU)."
        case default
            !print *, "OpenACC: Unknown or unsupported device type."
        end select
        
        !print *, "Starting molecule initialization"
        call init_molecule(mol, natoms, nelectrons, charge, nalpha, nbeta, z_num, coords)
        !print *, "Molecule initialized successfully."

        ! Call basis creation
        call build_basis(species, z_num, coords, natoms, basis_file, mol, shells)

        !call print_molecule(mol, shells)
        !print *, "Basis functions filled successfully."

        ! Run SCF calculation
        !if (nalpha == nbeta) then
        call run_rscf(mol, shells, eps_scf, max_iter, final_energy)
        !else 
            !call run_uscf(mol, shells, eps_scf, max_iter, final_energy)
        !end if
        !print *, "SCF calculation completed successfully."
        
    end subroutine run_scf_interface_py

end module main
