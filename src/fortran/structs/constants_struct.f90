module constants_struct
    use iso_c_binding
    implicit none

    integer, parameter :: wp = selected_real_kind(15)
    real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
    real(wp), parameter :: tpi = 2.0_wp*pi
    real(wp), parameter :: twopi25 = 2.0_wp*pi**(2.5_wp)
    real(wp), parameter :: top = 0.5_wp/atan(1.0_wp)

end module constants_struct