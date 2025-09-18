module molecule
    use iso_c_binding
    use basis_struct
    use openacc
    implicit none

contains

    subroutine init_molecule(mol, natoms, nelectrons, charge, nalpha, nbeta, atomic_numbers, coords)
        use iso_c_binding
        implicit none

        type(Molecule_type), intent(out) :: mol
        integer(c_int), intent(in) :: natoms, nalpha, nbeta, nelectrons, charge
        integer(c_int), intent(in) :: atomic_numbers(:)
        real(c_double), dimension(3, natoms), intent(in) :: coords        ! Expecting (3, natoms)


        integer :: ij, n, i, j, count
        real(wp) :: Rij, enuc_local

        ! Initialize scalar properties
        mol%n_atoms = natoms
        mol%nalpha = nalpha
        mol%nbeta = nbeta
        mol%nelectrons = nelectrons
        mol%charge = charge


        ! Initialize nuclear repulsion energy
        enuc_local = 0.0_wp

        ! OpenACC parallelization with proper array handling
        n = natoms * (natoms - 1) / 2
        !$acc parallel loop reduction(+:enuc_local)
        do ij = 1, n
            count = 0
            do i = 1, natoms - 1
                do j = i + 1, natoms
                    count = count + 1
                    if (count == ij) then
                        Rij = sqrt(sum((coords(:, i) - coords(:,j))**2))
                        if (Rij > 1.0e-10_wp) then
                            enuc_local = enuc_local + atomic_numbers(i) * atomic_numbers(j) / Rij
                        end if
                    end if
                end do
            end do
        end do

        ! Store result
        mol%Enuc = enuc_local

    end subroutine init_molecule

end module molecule
