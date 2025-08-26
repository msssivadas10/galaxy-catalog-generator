module hmf_tinker08_mod
    !! Halo mass-function model and Tinker (2008) and related models.

    use iso_c_binding
    ! use hmf_baselib_mod
    implicit none
    
    ! private

    integer, parameter :: dp = c_double

    real(c_double), parameter :: t08_table(5, 9) = reshape([ &
    !   Delta   ,    A       ,    a       ,    b       ,    c
         200._dp,    0.186_dp,    1.470_dp,    2.570_dp,    1.190_dp, & 
         300._dp,    0.200_dp,    1.520_dp,    2.250_dp,    1.270_dp, & 
         400._dp,    0.212_dp,    1.560_dp,    2.050_dp,    1.340_dp, & 
         600._dp,    0.218_dp,    1.610_dp,    1.870_dp,    1.450_dp, &
         800._dp,    0.248_dp,    1.870_dp,    1.590_dp,    1.580_dp, &
        1200._dp,    0.255_dp,    2.130_dp,    1.510_dp,    1.800_dp, &
        1600._dp,    0.260_dp,    2.300_dp,    1.460_dp,    1.970_dp, &
        2400._dp,    0.260_dp,    2.530_dp,    1.440_dp,    2.240_dp, &
        3200._dp,    0.260_dp,    2.660_dp,    1.410_dp,    2.440_dp  &
    ], shape(t08_table))  
    !! Table of Tinker (2008) mass-function parameters as function of overdensity Delta    

contains
    
end module hmf_tinker08_mod

program main
    use iso_c_binding
    use hmf_tinker08_mod
    implicit none

    integer :: i
    do i = 1, 9
        print *, t08_table(:, i)
    end do
    
end program main