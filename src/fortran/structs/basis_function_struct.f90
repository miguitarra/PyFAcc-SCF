module basis_function_struct
    use iso_c_binding
    use molecule_module
    !> Always declare everything explicitly
    implicit none

    type :: BasisFunction
        integer(c_int) :: atom
        real(c_double), allocatable :: exponents(:)
        real(c_double), allocatable :: coefficients(:)
        real(c_double), allocatable :: center(:) ! (3,)
    end type BasisFunction

end module basis_function_struct
