module main
    use iso_c_binding
    use molecule_struct
    use basis_function_struct
    use basis
    use constants_struct
    use molecule
    use openacc
    use rscf
    use uscf
    implicit none

contains

    subroutine run_scf_interface_py(natoms, nelectrons, nalpha, nbeta, species, coords, eps_scf, max_iter, final_energy, basis_set_n, z_num) bind(C, name="run_scf_interface_py")

        integer(c_int32_t), intent(in) :: natoms, nelectrons, nalpha, nbeta, max_iter, basis_set_n
        character(len=2), dimension(natoms), intent(in) :: species
        real(c_double), intent(in) :: coords(natoms, 3)
        integer(c_int32_t), intent(in) :: z_num(natoms)   
        real(c_double), intent(in) :: eps_scf
        real(c_double), intent(out) :: final_energy
        type(BasisFunction_type), allocatable :: basis_functions(:)
        type(Molecule_type) :: mol
        type(Shell_type), allocatable :: shells(:)
        integer :: dev_type
        integer :: implicit

        character(len=100) :: basis_file
        write(basis_file, '(A,I0,A)') './data/basis_set/sto-', basis_set_n, 'g.txt'

            
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
        call init_molecule(mol, natoms, nelectrons, nalpha, nbeta, z_num, coords)
        !print *, "Molecule initialized successfully."

        ! Call basis creation
        call build_basis_functions(species, coords, natoms, basis_file, basis_functions, shells)
        !print *, "Basis functions filled successfully."

        ! Run SCF calculation
        if (nalpha == nbeta) then
            call run_rscf(mol, shells, basis_functions, eps_scf, max_iter, final_energy, basis_set_n)
        else 
            call run_uscf(mol, shells, basis_functions, eps_scf, max_iter, final_energy, basis_set_n)
        end if
        !print *, "SCF calculation completed successfully."
        

    end subroutine run_scf_interface_py

end module main
