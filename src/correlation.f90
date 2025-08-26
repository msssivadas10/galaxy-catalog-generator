module correlation_mod
    !! Calculation of matter correlation function from power spectrum data.

    use, intrinsic :: ieee_arithmetic
    use iso_c_binding
    implicit none

    private
    public :: correlation

    real(c_double), parameter :: pi  = 3.141592653589793_c_double

contains

    function correlation_integrand(lnr, lnk, lnp, j) result(res)
        !! Calculate the value of integrand for variance integral.

        real(c_double), intent(in), value :: lnr
        !! Natural log of smoothing radius in Mpc
        
        real(c_double), intent(in) :: lnk
        !! Natural log of wavenumber in 1/Mpc

        real(c_double), intent(in) :: lnp
        !! Natural log of power spectrum value

        integer(c_int), intent(in), value :: j
        !! Non-zero value for j will give average matter correlation

        real(c_double) :: res, kr

        ! Dimenssionless matter power spectrum
        res = exp( lnk * (3 + 2*j) + lnp ) 
        
        ! Filter function calculation
        kr  = exp( lnk + lnr )
        if ( j == 0 ) then 
            res = res * ( sin(kr) / kr )
        else
            ! For average correlation function: this integral converges faster
            res = res * ( 3*( sin(kr) - kr * cos(kr) ) / kr**3 )
        end if

    end function correlation_integrand

    function correlation(lnr, j, pktab, size, cls) result(res) bind(C)
        !! Calculate the value of matter correlation function.

        real(c_double), intent(in), value :: lnr
        !! Natural log of smoothing radius in Mpc

        integer(c_int), intent(in), value :: j
        !! Non-zero value for j will give average matter correlation
        
        real(c_double), intent(in) :: pktab(cls, size)
        !! Precalculated power spectrum table. The nodes (lnk, first column) 
        !! should correspond to that of 7-th order Gauss-Legendre quadrature. 
        !! Main interval can be subdivided further and the nodes in each sub-
        !! interval can be stacked together, for more accuracy. 

        integer(c_int64_t), intent(in), value :: size
        !! Size of the power spectrum table: must be a multiple of 7.

        integer(c_int), intent(in), value :: cls
        !! Columns of the power spectrum table: must be 2 or 3. If 3, use the 
        !! the last column as weights. Otherwise, Simpson's rule is used. For 
        !! other values, return value will be NaN.

        real(c_double) :: res
        !! Calculated spectral moment value  

        integer(c_int64_t) :: i, m
        real(c_double)     :: lnk, fk, s, delta(2), weight(3)

        res = 0.0_c_double
        if ( cls == 2 ) then
            !! Use Simpson's rule integration
            
            do i = 1, size-2, 2
                ! For odd number of intervals, this loops runs up to the second
                ! last interval... 

                ! Calculating weights:
                delta = pktab(1, i+1:i+2) - pktab(1, i:i+1)
                if( abs( delta(2) - delta(1) ) < 1.0e-08_c_double ) then
                    ! Weights for uniform Simpson's rule
                    weight = [ 1.0_c_double, 4.0_c_double, 1.0_c_double ]
                else
                    ! Weights for non-uniform Simpson's rule
                    weight = [ &
                        2.0_c_double - delta(2) / delta(1), &
                        sum(delta)**2 / product(delta),     &
                        2.0_c_double - delta(1) / delta(2)  &
                    ]
                end if

                s = 0.0_c_double
                do m = i, i+2
                    lnk = pktab(1, m)
                    fk  = pktab(2, m)
                    fk  = correlation_integrand(lnr, lnk, fk, j) 
                    s   = s + fk * weight(m-i+1)
                end do
                res = res + sum(delta) * s
            
            end do 
            ! For odd number of intervals, the last two intervals are handled 
            ! seperately (ref:  <https://en.wikipedia.org/wiki/Simpson%27s_rule>)
            if ( modulo(size-1, 2) == 1 ) then
                i = size-2

                ! Calculating weights:
                delta = pktab(1, i+1:i+2) - pktab(1, i:i+1)
                if( abs( delta(2) - delta(1) ) < 1.0e-08_c_double ) then
                    ! Weights for uniform Simpson's rule
                    weight = [ 1.0_c_double, 4.0_c_double, 1.0_c_double ]
                else
                    ! Weights for non-uniform Simpson's rule
                    s      = sum(delta)
                    weight = [ &
                        delta(2) * ( 2*delta(2) + 3*delta(1) ) / sum(delta), &
                        delta(2) * (   delta(2) + 3*delta(1) ) / delta(1),   &
                        delta(2)**3 /  delta(1) / sum(delta)                 &
                    ]
                end if

                s = 0.0_c_double
                do m = i, i+2
                    lnk = pktab(1, m)
                    fk  = pktab(2, m)
                    fk  = correlation_integrand(lnr, lnk, fk, j) 
                    s   = s + fk * weight(m-i+1)
                end do
                res = res + sum(delta) * s

            end if
            res = res / 6.0_c_double

        else if ( cls == 3 ) then
            !! Use integration with col=3 as the weights
            do i = 1, size
                lnk = pktab(1, i)
                fk  = pktab(2, i)
                fk  = correlation_integrand(lnr, lnk, fk, j) 
                res = res + fk * pktab(3, i)
            end do

        else
            !! Invalid value: result will be NaN
            res = ieee_value(res, ieee_quiet_nan)
            return
        end if
        res = res / ( 2*pi**2 )
        
    end function correlation
    
end module correlation_mod