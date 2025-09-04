module halo_model_mod
    !! A module for basic halo model calculations. These calculations are based  
    !! on a 5-parameter halo occupation distribution (HOD) model. 

    use iso_fortran_env, only: stderr => error_unit
    use iso_c_binding
    use constants_mod
    use random_mod
    use interpolate_mod
    use integrate_mod
    implicit none

    private
    public :: lagrangian_r, central_count, satellite_count, subhalo_mass_function,  &
              halo_concentration, average_halo_density, average_galaxy_density,     &
              average_satellite_frac, average_galaxy_bias

    type, public, bind(c) :: hmargs_t
        !! A struct storing various halomodel parameters
        real(c_double) :: lnm_min    !! Minimum halo mass (in Msun) to have at least one central galaxy 
        real(c_double) :: sigma_m    !! Width of the central galaxy transition range. (0 for a step function)
        real(c_double) :: lnm0       !! Minimum halo mass (in Msun) to have satellite galaxies
        real(c_double) :: lnm1       !! Scale factor for power law satelliete count relation (Msun)
        real(c_double) :: alpha      !! Index for the  power law satelliete count relation
        real(c_double) :: scale_shmf !! Scale parameter for the subhalo mass-function
        real(c_double) :: slope_shmf !! Slope parameter for the subhalo mass-function
        real(c_double) :: z          !! Redshift
        real(c_double) :: H0         !! Hubble parameter value
        real(c_double) :: Om0        !! Total matter density parameter
        real(c_double) :: Delta_m    !! Matter overdensity w.r.to mean background density
        real(c_double) :: dplus      !! Growth factor at this redshift
    end type

contains

    function lagrangian_r(args, lnm) result(lnr) bind(c)
        !! Return the Lagrangian radius a halo (natural log of value in Mpc), 
        !! given its mass. 

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnm
        !! Natural log of halo mass (Msun)

        real(c_double) :: lnr, rho_m, rho_h

        rho_m = args%Om0 * ( critical_density_const * args%H0**2 ) ! Matter density at z=0 in Msun/Mpc^3 
        rho_h = rho_m ! Halo density (TODO: chek if the halo density is rho_m * self.Delta)

        ! Lagrangian radius (r) corresponding to halo mass
        lnr = ( lnm + log(3._c_double / (4*pi) / rho_h ) ) / 3._c_double ! r in Mpc
        
    end function lagrangian_r

    function central_count(args, lnm) result(res) bind(c)
        !! Return the average count of central galaxies in a halo, given its 
        !! mass. This will be a sigmoid function with smoothness controlled by 
        !! the `sigma_m` parameter. If it is 0, then it will be a step function.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnm
        !! Natural log of halo mass (Msun)

        real(c_double) :: res
    
        res = lnm - args%lnm_min

        if ( abs(args%sigma_m) < 1e-06 ) then
            ! Heaviside step function
            if ( res < 0._c_double ) then 
                res = 0._c_double
            else
                res = 1._c_double
            end if
        else
            ! Sigmoid function
            res = 0.5_c_double*( 1._c_double + erf(res / args%sigma_m) )
        end if

    end function central_count

    function satellite_count(args, lnm) result(res) bind(c)
        !! Return the average count of satellite galaxies in a halo, given 
        !! its mass. Average fraction of satellites is given by a power law 
        !! of the form `((m0 - m)/ m1)^alpha`.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnm
        !! Natural log of halo mass (Msun)

        real(c_double) :: res
        
        res = ( exp(lnm) - exp(args%lnm0) ) / exp(args%lnm1)
        if ( res < 0._c_double ) then
            res = 0._c_double
            return
        else
            res = res**args%alpha * central_count(args, lnm)
        end if     

        end function satellite_count

    function subhalo_mass_function(args, x, lnm) result(res) bind(c)
        !! Calculate the subhalo mass-function for given halo mass. This 
        !! is a bounded power law defined in the subhalo masses in the 
        !! range `[m_min, scale*m]`.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: x
        !! Mass of the subhalo as fraction of the parent halo

        real(c_double), intent(in), value :: lnm
        !! Natural log of halo mass (Msun)

        real(c_double) :: res, a, b, c, p
        
        a = exp(args%lnm_min - lnm)
        b = args%scale_shmf
        if ( x < a .or. x > b ) then
            res = 0._c_double
            return
        end if
        
        p   = args%slope_shmf ! power law index
        c   = p * a**p / (1._c_double - (a / b)**p) ! amplitude
        res = c * x**(-p - 1)
    
    end function subhalo_mass_function

    function halo_concentration(args, sigma) result(res) bind(c)
        !! Return the value of halo concentration parameter for a given 
        !! mass, calculated for the current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: sigma
        !! Matter variance corresponding to the halo mass

        real(c_double) :: res, zp1, c0, b, g1, g2, v0, v, t

        ! Redshift dependent parameters:
        zp1 = (1._c_double + args%z)
        c0  = 3.395_c_double * zp1**(-0.215_c_double )
        b   = 0.307_c_double * zp1**( 0.540_c_double )
        g1  = 0.628_c_double * zp1**(-0.047_c_double )
        g2  = 0.317_c_double * zp1**(-0.893_c_double )
        v0  = ( &
                4.135_c_double                   &
                    - 0.564_c_double   * zp1     &
                    - 0.210_c_double   * zp1**2  &
                    + 0.0557_c_double  * zp1**3  &
                    - 0.00348_c_double * zp1**4  &
             ) / args%dplus
        v  = delta_sc / sigma
        t  = v / v0

        ! Concentration-mass relation:
        res = c0 * t**(-g1) * (1._c_double + t**(1._c_double / b))**(-b*(g2 - g1))
        
    end function halo_concentration

! Halo averages:

    function average_halo_density(args, lnma, lnmb, mfspline, mfns, &
                                  abstol, reltol, maxiter           &
        ) result(res) bind(c)
        !! Return the average halo number density for the given halo mass 
        !! range at current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnma
        !! Lower limit: natural log of halo mass (Msun)
        
        real(c_double), intent(in), value :: lnmb
        !! Upper limit: natural log of halo mass (Msun)

        real(c_double), intent(in) :: mfspline(3,mfns)
        !! Cubic spline for halo mass function, log(dn/dlnm) as function 
        !! of log(halo mass)

        integer(c_int64_t), intent(in), value :: mfns
        !! Size of the mass-function spline data

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double) :: res, err, empty(3, 0) 
        integer(c_int) :: stat
        
        call calc_halo_average(0, lnma, lnmb, args, mfspline, mfns, empty, 0_c_int64_t, &
                               abstol, reltol, maxiter, res, err, stat                  &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_halo_density: integral failed to converge'
        
    end function average_halo_density

    function average_galaxy_density(args, lnma, lnmb, mfspline, mfns, &
                                    abstol, reltol, maxiter           &
        ) result(res) bind(c)
        !! Return the average galaxy number density for the given halo mass 
        !! range at current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnma
        !! Lower limit: natural log of halo mass (Msun)
        
        real(c_double), intent(in), value :: lnmb
        !! Upper limit: natural log of halo mass (Msun)

        real(c_double), intent(in) :: mfspline(3,mfns)
        !! Cubic spline for halo mass function, log(dn/dlnm) as function 
        !! of log(halo mass)

        integer(c_int64_t), intent(in), value :: mfns
        !! Size of the mass-function spline data

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double) :: res, res1, res2, err, empty(3, 0) 
        integer(c_int) :: stat
        
        ! Average central galaxy density
        call calc_halo_average(1, lnma, lnmb, args, mfspline, mfns, empty, 0_c_int64_t, &
                               abstol, reltol, maxiter, res1, err, stat                 &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_galaxy_density: integral (c) failed to converge'
        
        ! Average satellite galaxy density
        call calc_halo_average(2, lnma, lnmb, args, mfspline, mfns, empty, 0_c_int64_t, &
                               abstol, reltol, maxiter, res2, err, stat                 &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_galaxy_density: integral (s) failed to converge'

        ! Average galaxy density
        res = res1 + res2

    end function average_galaxy_density

    function average_satellite_frac(args, lnma, lnmb, mfspline, mfns, &
                                    abstol, reltol, maxiter           &
        ) result(res) bind(c)
        !! Return the average satellite fraction for the given halo mass 
        !! range at current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnma
        !! Lower limit: natural log of halo mass (Msun)
        
        real(c_double), intent(in), value :: lnmb
        !! Upper limit: natural log of halo mass (Msun)

        real(c_double), intent(in) :: mfspline(3,mfns)
        !! Cubic spline for halo mass function, log(dn/dlnm) as function 
        !! of log(halo mass)

        integer(c_int64_t), intent(in), value :: mfns
        !! Size of the mass-function spline data

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double) :: res, res1, res2, err, empty(3, 0) 
        integer(c_int) :: stat
        
        ! Average central galaxy density
        call calc_halo_average(1, lnma, lnmb, args, mfspline, mfns, empty, 0_c_int64_t, &
                               abstol, reltol, maxiter, res1, err, stat                 &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_satellite_frac: integral (c) failed to converge'
        
        ! Average satellite galaxy density
        call calc_halo_average(2, lnma, lnmb, args, mfspline, mfns, empty, 0_c_int64_t, &
                               abstol, reltol, maxiter, res2, err, stat                 &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_satellite_frac: integral (s) failed to converge'

        ! Average satellite fraction
        res = res2 / (res1 + res2)
        
    end function average_satellite_frac

    function average_galaxy_bias(args, lnma, lnmb, mfspline, mfns, bfspline, &
                                 bfns, abstol, reltol, maxiter               &
        ) result(res) bind(c)
        !! Return the average satellite fraction for the given halo mass 
        !! range at current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnma
        !! Lower limit: natural log of halo mass (Msun)

        real(c_double), intent(in), value :: lnmb
        !! Upper limit: natural log of halo mass (Msun)

        real(c_double), intent(in) :: mfspline(3,mfns)
        !! Cubic spline for halo mass function, log(dn/dlnm) as function 
        !! of log(halo mass)

        integer(c_int64_t), intent(in), value :: mfns
        !! Size of the mass-function spline data

        real(c_double), intent(in) :: bfspline(3,bfns)
        !! Bias spline

        integer(c_int64_t), intent(in), value :: bfns
        !! Size of the bias spline

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double) :: res, res1, res2, res3, res4, err, empty(3, 0)
        integer(c_int) :: stat

        ! Average central galaxy density
        call calc_halo_average(1, lnma, lnmb, args, mfspline, mfns, empty, 0_c_int64_t, &
                               abstol, reltol, maxiter, res1, err, stat                 &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_galaxy_bias: integral (c) failed to converge'

        ! Average satellite galaxy density
        call calc_halo_average(2, lnma, lnmb, args, mfspline, mfns, empty, 0_c_int64_t, &
                               abstol, reltol, maxiter, res2, err, stat                 &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_galaxy_bias: integral (s) failed to converge'

        ! Average central galaxy bias
        call calc_halo_average(1, lnma, lnmb, args, mfspline, mfns, bfspline, bfns, &
                               abstol, reltol, maxiter, res3, err, stat             &
        )
        if ( stat /= 0 ) &
        write(stderr,'(a)') 'warning: average_galaxy_bias: integral (bc) failed to converge'

        ! Average satellite galaxy bias
        call calc_halo_average(2, lnma, lnmb, args, mfspline, mfns, bfspline, bfns, &
                               abstol, reltol, maxiter, res4, err, stat             &
        )
        if ( stat /= 0 ) &
        write(stderr,'(a)') 'warning: average_galaxy_bias: integral (bs) failed to converge'

        ! Average galaxy bias
        res = (res3 + res4) / (res1 + res2)

        end function average_galaxy_bias

    subroutine calc_halo_average(f, a, b, args, mfspline, mfns, bfspline, bfns, &
                                 abstol, reltol, maxiter, res, err, stat        &
        )
        !! Calculate the halo average with optional bias weight. 

        integer(c_int), intent(in), value :: f
        !! Function selector (0=no weight, 1=central count, 2=satellite count)

        real(c_double), intent(in), value :: a
        !! Lower limit of integration

        real(c_double), intent(in), value :: b
        !! Upper limit of integration

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in) :: mfspline(3,mfns)
        !! Mass-function spline

        integer(c_int64_t), intent(in), value :: mfns
        !! Size of the mass-function spline

        real(c_double), intent(in) :: bfspline(3,bfns)
        !! Bias spline

        integer(c_int64_t), intent(in), value :: bfns
        !! Size of the bias spline

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double), intent(out) :: res
        !! Value of the integral of f over [a, b]

        real(c_double), intent(out) :: err
        !! Estimate of the error in integration

        integer(c_int), intent(out) :: stat
        !! Error code: 0=ok, 1=integral not converged

        integer(c_int64_t) :: iter
        real(c_double)     :: xa, xb, xm, I0, I1, I2, err0, err1, err2
        integer(c_int64_t) :: heap_size, heap_capacity
        real(c_double), allocatable :: heap(:, :)

        heap_size     = 0
        heap_capacity = 10*maxiter ! Heap capacity
        allocate( heap(4, heap_capacity) )

        ! Initial evaluation
        call calc_halo_average2(f, a, b, args, mfspline, mfns, &
                                bfspline, bfns, I0, err0       &
        )
        call int_heap_push(heap, heap_size, heap_capacity, a, b, I0, err0)

        res  = I0
        err  = err0
        stat = 1
        do iter = 1, maxiter

            ! Stop if tolerance is met
            if ( err <= max(abstol, reltol*abs(res)) ) then
                stat = 0 
                exit
            endif

            ! Pop worst interval
            call int_heap_pop(heap, heap_size, heap_capacity, xa, xb, I0, err0)

            xm = 0.5_c_double * (xa + xb)
            
            ! Refine on left interval
            call calc_halo_average2(f, xa, xm, args, mfspline, mfns, &
                                    bfspline, bfns, I1, err1         &
            )
            call int_heap_push(heap, heap_size, heap_capacity, xa, xm, I1, err1) ! Push new interval back
            
            ! Refine on left interval
            call calc_halo_average2(f, xm, xb, args, mfspline, mfns, &
                                    bfspline, bfns, I2, err2         &
            )
            call int_heap_push(heap, heap_size, heap_capacity, xm, xb, I2, err2) ! Push new interval back
            
            ! Update global sums
            res = res + (I1   + I2   - I0  ) ! replace old interval
            err = err + (err1 + err2 - err0)
            
        end do

        deallocate(heap)    
        
    end subroutine calc_halo_average

    subroutine calc_halo_average2(f, a, b, args, mfspline, mfns, bfspline, &
                                  bfns, res, err                           &
        )
        !! Calculate the halo average with optional bias weight in the 
        !! interval [a, b]. 

        integer(c_int), intent(in), value :: f
        !! Function selector (0=no weight, 1=central count, 2=satellite count)

        real(c_double), intent(in), value :: a
        !! Lower limit of integration

        real(c_double), intent(in), value :: b
        !! Upper limit of integration

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in) :: mfspline(3,mfns)
        !! Mass-function spline

        integer(c_int64_t), intent(in), value :: mfns
        !! Size of the mass-function spline

        real(c_double), intent(in) :: bfspline(3,bfns)
        !! Bias spline

        integer(c_int64_t), intent(in), value :: bfns
        !! Size of the bias spline

        real(c_double), intent(out) :: res
        !! Value of the integral of f over [a, b]

        real(c_double), intent(out) :: err
        !! Estimate of the error in integration
            
        integer(c_int64_t) :: j
        real(c_double)     :: intg, intk, xval, fval, fval2, scale

        scale = 0.5_c_double * (b - a)
        xval  = a + scale ! log(m)
        fval  = exp( interpolate(xval, mfns, mfspline) ) ! dn/dlnm
        if ( bfns > 0 ) fval = fval * exp( interpolate(xval, bfns, bfspline) ) ! bias
        select case ( f )
        case ( 1 ) ! weight=central galaxy count
            fval = fval * central_count(args, xval)
        case ( 2 ) ! weight=satellite galaxy count
            fval = fval * satellite_count(args, xval)
        ! default=no weight
        end select 
        intk  = fval * K15(2,1) 
        intg  = fval *  G7(2,1) 
        do j = 2, 8

            xval = a + scale * (1. - K15(1,j)) ! log(m)
            fval = exp( interpolate(xval, mfns, mfspline) ) ! dn/dlnm
            if ( bfns > 0 ) fval = fval * exp( interpolate(xval, bfns, bfspline) ) ! bias
            select case ( f )
            case ( 1 ) ! weight=central galaxy count
                fval = fval * central_count(args, xval)
            case ( 2 ) ! weight=satellite galaxy count
                fval = fval * satellite_count(args, xval)
            ! default=no weight
            end select 

            xval  = a + scale * (1. + K15(1,j)) ! log(m)
            fval2 = exp( interpolate(xval, mfns, mfspline) ) ! dn/dlnm
            if ( bfns > 0 ) fval2 = fval2 * exp( interpolate(xval, bfns, bfspline) ) ! bias
            select case ( f )
            case ( 1 ) ! weight=central galaxy count
                fval2 = fval2 * central_count(args, xval)
            case ( 2 ) ! weight=satellite galaxy count
                fval2 = fval2 * satellite_count(args, xval)
            ! default=no weight
            end select 
            fval = fval + fval2
            
            intk = intk + fval * K15(2,j)
            if ( mod(j, 2) == 1 ) intg = intg + fval * G7(2,(j+1)/2) ! Point also in G7 rule
        end do
        intk = scale * intk
        intg = scale * intg
        res  = intk
        err  = abs(intk - intg)
        
    end subroutine calc_halo_average2

end module halo_model_mod
