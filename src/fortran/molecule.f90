module molecule
    use iso_c_binding
    use molecule_struct
    use openacc
    implicit none

contains

    subroutine init_molecule(mol, natoms, nelectrons, nalpha, nbeta, atomic_numbers, coords)
        use iso_c_binding
        use molecule_struct
        implicit none

        type(Molecule_type), intent(out) :: mol
        integer(c_int), intent(in) :: natoms, nalpha, nbeta, nelectrons
        integer(c_int), intent(in) :: atomic_numbers(:)
        real(c_double), dimension(3, natoms), intent(in) :: coords        ! Expecting (3, natoms)


        integer :: ij, n, i, j, count
        real(wp) :: Rij, enuc_local

        ! Initialize scalar properties
        mol%natoms = natoms
        mol%nalpha = nalpha
        mol%nbeta = nbeta
        mol%nelectrons = nelectrons

        ! Allocate and assign molecule data
        allocate(mol%atomic_numbers(natoms))
        allocate(mol%coords(3, natoms))

        mol%atomic_numbers = atomic_numbers
        mol%coords = coords


        ! Initialize nuclear repulsion energy
        enuc_local = 0.0_wp

        ! OpenACC parallelization with proper array handling
        n = mol%natoms * (mol%natoms - 1) / 2
        !$acc parallel loop reduction(+:enuc_local)
        do ij = 1, n
            count = 0
            do i = 1, mol%natoms - 1
                do j = i + 1, mol%natoms
                    count = count + 1
                    if (count == ij) then
                        Rij = sqrt(sum((mol%coords(:, i) - mol%coords(:,j))**2))
                        if (Rij > 1.0e-10_wp) then
                            enuc_local = enuc_local + mol%atomic_numbers(i) * mol%atomic_numbers(j) / Rij
                        end if
                    end if
                end do
            end do
        end do

        ! Store result
        mol%Enuc = enuc_local

    end subroutine init_molecule

end module molecule
