module hmf_tinker08_mod
    !! Halo mass-function model and Tinker (2008) and related models.

    use iso_c_binding
    use hmf_baselib_mod
    implicit none
    
    ! private

    integer, parameter :: dp = c_double

    real(c_double), parameter :: T08(5, 9) = reshape([ &
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
    ], shape(T08))  
    !! Table of Tinker (2008) mass-function parameters as function of overdensity Delta    

contains

    subroutine setup_hmf_tinker08(args) bind(c)
        !! Setup for Tinker(2008) halo mass-function calculation.
        
        type(hmfargs_t), intent(inout) :: args

        ! real(c_double) :: d(4,9), u(4,8), s, dx1, dx2, dy1(4), dy2(4), p(4)
        integer :: k, klo, khi

        ! Get the values for parameters: interpolate from the table

        !! -- Cubic interpolation
        ! d(:,1) = 0.0_c_double
        ! d(:,9) = 0.0_c_double
        ! u(:,1) = 0.0_c_double
        ! do k = 2, 8
        !     dx1    = T08(1  , k  ) - T08(1  , k-1)
        !     dx2    = T08(1  , k+1) - T08(1  , k  )
        !     dy1    = T08(2:5, k  ) - T08(2:5, k-1)
        !     dy2    = T08(2:5, k+1) - T08(2:5, k  )
        !     s      = dx1 / (dx1 + dx2)
        !     p      = s * d(:,k-1) + 2._c_double
        !     u(:,k) = 6*( dy2/dx2 - dy1/dx1 ) / (dx1 + dx2) - s*u(:,k-1)
        ! end do
        klo = 1
        khi = 9
        do while( khi - klo > 1 )
            k = (khi + klo) / 2
            if ( T08(1,k) > args%Delta_m ) then
                khi = k
            else
                klo = k
            end if
        end do
        ! s   = T08(1,khi) - T08(1,klo)
        ! dx1 = ( T08(1,khi) - args%Delta_m     ) / s
        ! dx1 = ( args%Delta_m     - T08(1,klo) ) / s
        ! p   = dx1*T08(2:5,klo) + dx2*T08(2:5,khi)
        ! dx1 = ( dx1**3 - dx1 )
        ! dx2 = ( dx2**3 - dx2 )
        ! p   = p + (dx1 * d(:,klo) + dx2 * d(:,khi)) * s**2 / 6._c_double
        ! args%param(1:4) = p(:) ! Parameters [ A0, a0, b0, c0 ]

        ! -- Linear interpolation:
        args%param(1:4) = (T08(2:5,khi) - T08(2:5,klo)) / (T08(1,khi) - T08(1,klo)) ! slope
        args%param(1:4) = T08(2:5,klo) + args%param(1:4) * (args%Delta_m - T08(1,klo))
        ! Parameters [ A0, a0, b0, c0 ]

        args%param(5) = -( 0.75_c_double / log10( args%Delta_m / 75._c_double ) )**1.2_c_double 
        args%param(5) =  10._c_double**args%param(5)
        
    end subroutine setup_hmf_tinker08

    function hmf_tinker08(args, target_code) result(res)
        !! Calculate Tinker (2008) halo mass-unssfunction value
        
        type(hmfargs_t), intent(in) :: args
        !! Arguments

        integer(c_int), intent(in), value :: target_code
        !! Code for output value (0=dn/dm, 1=dn/dlnm, 2=dn/dlog10m)
        
        real(c_double) :: res, A_, a, b, c, alpha, zp1, s

        zp1 = args%z + 1
        
        ! Parameter `A`
        A_ = args%param(1) * zp1**(-0.14_c_double ) 
        
        ! Parameter `a`  
        a = args%param(2) * zp1**(-0.06_c_double ) 
        
        ! Parameter `b`
        alpha = args%param(5)
        b     = args%param(3) * zp1**(-alpha)
        
        ! Parameter `c`
        c = args%param(4)
        
        ! Tinker (2008) mass-function, f(s)
        s   = args%s
        res = A_ * (1._c_double + (b / s)**a) * exp(-c / s**2)
        
    end function hmf_tinker08
    
end module hmf_tinker08_mod