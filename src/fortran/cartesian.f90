module cartesian
    use basis_struct
    implicit none

contains

    subroutine generate_cartesian_angmoms(l, angmoms)
        !$acc routine seq
        integer, intent(in) :: l
        integer, intent(out) :: angmoms(:,:)
        integer :: ax, ay, az, idx, n
    
        idx = 1
        do ax = l, 0, -1
            do ay = l-ax, 0, -1
                az = l - ax - ay
                angmoms(idx,1) = ax
                angmoms(idx,2) = ay
                angmoms(idx,3) = az
                idx = idx + 1
            end do
        end do
    end subroutine generate_cartesian_angmoms

    integer function count_nbf(shells) result(nbf)
        ! input
        type(Shell_type), intent(in) :: shells(:)
    
        ! local
        integer :: nbf, n_shells, i
    
        n_shells = size(shells)
        nbf = 0
        do i = 1, n_shells
            nbf = nbf + ( (shells(i)%l_num + 1)*(shells(i)%l_num + 2)/2 )
        end do
        
    end function count_nbf

    

end module