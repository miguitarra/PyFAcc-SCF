module la
    use iso_c_binding
    use constants_struct
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
            !stop
        end if
    end subroutine diagonalize_symmetric

    subroutine orthonormal_diagonalize(F, S, eigvals, eigvecs)
        real(wp), intent(in)  :: F(:,:), S(:,:)
        real(wp), intent(out) :: eigvals(:), eigvecs(:,:)
        real(wp), allocatable :: S_copy(:,:), S_eigvals(:), S_eigvecs(:,:), X(:,:), F_prime(:,:)
        integer :: n, i, j

        n = size(F,1)
        allocate(S_copy(n,n), S_eigvals(n), S_eigvecs(n,n), X(n,n), F_prime(n,n))

        S_copy = S
        call diagonalize_symmetric(S_copy, S_eigvals, S_eigvecs)

        do j = 1, n
            if (S_eigvals(j) < 1.0e-8_wp) then
                print *, "Warning: Small overlap eigenvalue detected: ", S_eigvals(j)
            end if
        end do

        do i = 1, n
            do j = 1, n
                X(i,j) = S_eigvecs(i,j) / sqrt(S_eigvals(j))
            end do
        end do

        F_prime = matmul(transpose(X), matmul(F, X))
        call diagonalize_symmetric(F_prime, eigvals, eigvecs)
        eigvecs = matmul(X, eigvecs)

        deallocate(S_copy, S_eigvals, S_eigvecs, X, F_prime)
    end subroutine orthonormal_diagonalize


    subroutine idx_upper_triangular(n, idx, i, j)
        !$acc routine seq
        integer, intent(in) :: n, idx ! n for n x n matrix, idx linear index
        integer, intent(out) :: i, j ! i,j location
        integer :: k
    
        k = n
        i = 1
        j = idx
        do while (j > k)
            j = j - k
            k = k - 1
            i = i + 1
        end do
        j = i + j - 1
    end subroutine
end module la