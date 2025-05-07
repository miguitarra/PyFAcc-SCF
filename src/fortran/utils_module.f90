module utils_module
    use iso_c_binding
    use molecule_module
    implicit none
    
    contains
    
    subroutine diagonalize_symmetric(A, eigvals, eigvecs)
        real(c_double), intent(inout) :: A(:,:)
        real(c_double), intent(out) :: eigvals(:)
        real(c_double), intent(out) :: eigvecs(:,:)
    
        
        integer :: n, info
        real(c_double), allocatable :: work(:)

        n = size(A,1)
        allocate(work(3*n))
        A = 0.5_wp * (A + transpose(A))
        call dsyev('V', 'U', n, A, n, eigvals, work, 3*n, info)
        eigvecs = A

        if (info /= 0) then
            print *, "Diagonalization failed"
            stop
        end if
    end subroutine diagonalize_symmetric


    end module utils_module
     