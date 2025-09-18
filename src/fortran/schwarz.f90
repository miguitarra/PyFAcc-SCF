module schwarz
    use constants_struct
    use nuclear
    use utils
    use eri_struct
    use boys
    implicit none

    real(8), parameter :: threshold = 1.0D-10

contains 

    subroutine schwarz_matrix(mol, basis_function, Sch)
        !Input
        type(Molecule_type), intent(in) :: mol
        type(BasisFunction_type), intent(in) :: basis_function(:)

        


    end subroutine schwarz_matrix

end module schwarz