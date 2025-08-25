module power_integrals
    !!
    !! Functions for evaluating various integrals related to powerspectrum. 
    !!
    use iso_c_binding
    use quadutils
    implicit none

    private
    public :: variance, correlation

    real(c_double), parameter :: pi  = 3.141592653589793_c_double
    
contains

    function variance(lnr, j, nu, filt, pktab, size, rule) result(res) bind(C)
        !! Calculate the value of smoothed matter variance for j=0 and nu=0.
        !!
        !! Calculate the j-th order spectral moment of the density field, smoothed 
        !! using the specified filter function at a radius of r. 

        real(c_double), intent(in), value :: lnr
        !! Natural log of smoothing radius in Mpc

        integer(c_int), intent(in), value :: j
        !! Order of the moment (j=0 will give matter variance)

        integer(c_int), intent(in), value :: nu
        !! If non-zero, calculate first derivative w.r.to r

        integer(c_int), intent(in), value :: filt
        !! Code of filter function (0=tophat, 1=gaussian)
        
        real(c_double), intent(in) :: pktab(2, size)
        !! Precalculated power spectrum table. The nodes (lnk, first column) 
        !! should correspond to that of 7-th order Gauss-Legendre quadrature. 
        !! Main interval can be subdivided further and the nodes in each sub-
        !! interval can be stacked together, for more accuracy. 
        
        integer(c_int64_t), intent(in), value :: size
        !! Size of the power spectrum table: must be a multiple of 7.

        integer(c_int), intent(in), value :: rule
        !! Integration rule. 0=stacked G7, 1=Simpson's rule with non-uniform spacing,
        !! 2=Simpson'r rule with uniform spacing.

        real(c_double) :: res
        !! Calculated spectral moment value  

        integer(c_int) :: stat
        real(c_double), dimension(size) :: lnk, fk, kr, wk, dwk
        
        lnk = pktab(1, :)
        fk  = exp( pktab(2, :) + (3 + 2*j)*lnk ) ! Dimenssionless matter power spectrum
        
        ! Filter function calculation
        kr = exp(lnk + lnr)
        if ( filt == 1 ) then
            ! Gaussian filter function
            wk = exp(-0.5_c_double * kr**2)
        else
            ! Spherical tophat function
            wk = 3*( sin(kr) - kr * cos(kr) ) / kr**3
        end if
        
        if ( nu == 0 ) then
            fk = fk * wk * wk
        else
            ! If calculating derivative:
            if ( filt == 1 ) then
                ! Gaussian filter function
                dwk = -kr * wk
            else
                ! Spherical tophat function
                dwk = ( 3*(kr**2 - 3) * sin(kr) + 9*kr * cos(kr) ) / kr**4
            end if
            fk = fk * 2*wk * dwk * exp(lnk)
        end if
        
        call integrate_functable(lnk, fk, rule, size, res, stat)
        res = res / ( 2*pi**2 )
        
    end function variance

    function correlation(lnr, j, pktab, size, rule) result(res) bind(C)
        !! Calculate the value of matter correlation function.

        real(c_double), intent(in), value :: lnr
        !! Natural log of smoothing radius in Mpc

        integer(c_int), intent(in), value :: j
        !! Non-zero value for j will give average matter correlation
        
        real(c_double), intent(in) :: pktab(2, size)
        !! Precalculated power spectrum table. The nodes (lnk, first column) 
        !! should correspond to that of 7-th order Gauss-Legendre quadrature. 
        !! Main interval can be subdivided further and the nodes in each sub-
        !! interval can be stacked together, for more accuracy. 

        integer(c_int64_t), intent(in), value :: size
        !! Size of the power spectrum table: must be a multiple of 7.

        integer(c_int), intent(in), value :: rule
        !! Integration rule. 0=stacked G7, 1=Simpson's rule with non-uniform spacing,
        !! 2=Simpson'r rule with uniform spacing.

        real(c_double) :: res
        !! Calculated spectral moment value  

        integer(c_int) :: stat
        real(c_double), dimension(size) :: lnk, kr, fk
        
        lnk = pktab(1, :)
        fk  = exp( pktab(2, :) + 3*lnk ) ! Dimenssionless matter power spectrum

        ! Filter function calculation
        kr = exp(lnk + lnr)
        if ( j == 1 ) then
            ! For average correlation function: this integral converges faster
            fk = fk * ( 3*( sin(kr) - kr * cos(kr) ) / kr**3 )
        else
            fk = fk * sin(kr) / kr
        end if

        call integrate_functable(lnk, fk, rule, size, res, stat)
        res = res / ( 2*pi**2 )
        
    end function correlation
    
end module power_integrals
