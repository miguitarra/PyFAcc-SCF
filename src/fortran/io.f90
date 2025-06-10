module io
    use iso_c_binding
    use constants_struct
    implicit none
contains
    subroutine print_matrix(matrix, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: matrix(n, n)
        integer :: i, j
        
        do i = 1, n
        write(*,'(100(f12.6,1x))') (matrix(i, j), j = 1, n)
        end do
    end subroutine print_matrix

    subroutine print_vector(vec, n)
        implicit none
        integer, intent(in) :: n
        real(c_double), intent(in) :: vec(n)
        integer :: i
    
        write(*,'(100(f12.6,1x))') (vec(i), i = 1, n)
    end subroutine print_vector

    function itoa(i) result(str)
        integer, intent(in) :: i
        character(len=12) :: str
        
        write(str, '(I0)') i
    end function itoa

end module io