module utils
    use iso_c_binding
    use molecule_struct
    use openacc
    implicit none
    
contains
    
    pure function factorial2(n) result(fact2)
        !$acc routine seq

        ! INPUT
        integer, intent(in) :: n

        ! INTERMEDIATE VARIABLE
        integer :: i

        ! OUTPUT
        integer :: fact2

        if (n == -1) then
            fact2 = 1
        else if (MOD(n,2) == 0) then ! EVEN NUMBER
            fact2 = 1

            do i = 2, n, 2
                fact2 = fact2 * i
            end do
        else ! ODD NUMBER
            fact2 = 1

            do i = 1, n, 2
                fact2 = fact2 * i
            end do
        end if

    end function factorial2


    pure function factorial(n) result(res)
        !$acc routine seq
        integer, intent(in) :: n
        integer :: i

        ! OUTPUT
        integer :: res

        res = 1

        do i = 2, n
            res = res * i
        end do
    end function factorial


    function binom(n,k)
        !$acc routine seq

        ! INPUT
        integer, intent(in) :: n, k

        ! OUTPUT
        integer :: binom

        binom = factorial(n)
        binom = binom / factorial(k)
        binom = binom / factorial(n-k)

    end function binom


end module utils
     