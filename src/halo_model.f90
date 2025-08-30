module halo_model_mod
    !! Halo occupation distributions and related calculations.

    use iso_fortran_env, only: stderr => error_unit
    use iso_c_binding
    use random_mod
    use interpolate_mod
    use integrate_mod
    implicit none

    private
    public :: central_count, satellite_count, subhalo_mass_function,             &
              halo_concentration, setup_catalog_generation, generate_satellites, &
              average_halo_density, average_galaxy_density, average_satellite_frac

    real(c_double), parameter :: delta_sc = 1.6864701998411453_c_double
    !! Overdensity threshold for spherical collapse in EdS universe
    
    real(c_double), parameter :: critical_density_const = 2.775e+11_c_double
    !! Constant part of the present critical density
    
    real(c_double), parameter :: pi = 3.141592653589793_c_double
    !! Pi

    type, public, bind(c) :: hmargs_t
        !! A struct storing various halomodel parameters
        
        real(c_double) :: lnm_min
        !! Minimum halo mass (in Msun) to have at least one central galaxy 

        real(c_double) :: sigma_m
        !! Width of the central galaxy transition range. (0 for a step function)

        real(c_double) :: lnm0
        !! Minimum halo mass (in Msun) to have satellite galaxies
        
        real(c_double) :: lnm1
        !! Scale factor for power law satelliete count relation (Msun)

        real(c_double) :: alpha
        !! Index for the  power law satelliete count relation

        real(c_double) :: scale_shmf
        !! Scale parameter for the subhalo mass-function

        real(c_double) :: slope_shmf
        !! Slope parameter for the subhalo mass-function

        real(c_double) :: z
        !! Redshift
        
        real(c_double) :: H0
        !! Hubble parameter value
        
        real(c_double) :: Om0
        !! Total matter density parameter

        real(c_double) :: dplus
        !! Growth factor at this redshift

    end type

    type, public, bind(c) :: cgargs_t
        !! Struct containing various parameters for galaxy catalog
        !! generation (execpt halo model parameters). 

        real(c_double) :: lnm 
        !! Natural log of halo mass in Msun

        real(c_double) :: pos(3)
        !! Coordinates of the halo in Mpc units

        real(c_double) :: lnr
        !! Hatural log of halo radius in Mpc

        real(c_double) :: c
        !! Concentration parameter for the halo

        integer(c_int64_t) :: n_cen
        !! Number of central galaxies (0 or 1)

        integer(c_int64_t) :: n_sat
        !! Number of satellite galaxies

    end type cgargs_t
    
contains

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
        res = c * x**(p - 1)
    
    end function subhalo_mass_function

    function halo_concentration(args, lnm, s_spline, ns) result(res) bind(c)
        !! Return the value of halo concentration parameter for a given 
        !! mass, calculated for the current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnm
        !! Natural log of halo mass (Msun)

        real(c_double), intent(in) :: s_spline(3,ns)
        !! Cubic spline for matter variance, as function of log(halo mass)

        integer(c_int64_t), intent(in), value :: ns
        !! Size of the matter variance spline data

        real(c_double) :: res, a, c0, b, g1, g2, v0, v, t

        ! Redshift dependent parameters:
        a  = 1._c_double / (1._c_double + args%z)
        c0 = 3.395_c_double * a**( 0.215_c_double )
        b  = 0.307_c_double * a**(-0.540_c_double )
        g1 = 0.628_c_double * a**( 0.047_c_double )
        g2 = 0.317_c_double * a**( 0.893_c_double )
        v0 = ( &
                4.135_c_double                   &
                    - 0.564_c_double   / a       &
                    - 0.210_c_double   / a**2    &
                    + 0.0557_c_double  / a**3    &
                    - 0.00348_c_double / a**4    &
             ) / args%dplus
        v  = delta_sc / interpolate(lnm, ns, s_spline)
        t  = v / v0

        ! Concentration-mass relation:
        res = c0 * t**(-g1) * (1._c_double + t**(1._c_double / b))**(-b*(g2 - g1))
        
    end function halo_concentration

! Galaxy Catalog Generation: 
    
    subroutine setup_catalog_generation(params, s_spline, ns, args, rstate) bind(c)
        !! Calculate various parameters for a galaxy catalog generation

        type(hmargs_t), intent(in) :: params 
        !! Halo model parameters

        real(c_double), intent(in) :: s_spline(3,ns)
        !! Cubic spline for matter variance, as function of log(halo mass)

        integer(c_int64_t), intent(in), value :: ns
        !! Size of the matter variance spline data

        type(cgargs_t), intent(inout) :: args
        !! Arguments for catalog generation

        integer(c_int64_t), intent(inout) :: rstate
        !! Random number generator state

        real(c_double) :: rho_m, rho_h, p_cen, lam_sat

        rho_m = params%Om0 * ( critical_density_const * params%H0**2 ) ! Matter density at z=0 in Msun/Mpc^3 
        rho_h = rho_m ! Halo density (TODO: chek if the halo density is rho_m * self.Delta)

        ! Lagrangian radius (r) corresponding to halo mass
        args%lnr = ( args%lnm + log(3._c_double / (4*pi) / rho_h ) ) / 3._c_double ! r in Mpc

        ! Central galaxy count: this is drawn from a binomial distribution of n = 1, 
        ! so that the value is either 1 or 0. If the model uses a step function use 
        ! the average count as the actual count. 
        p_cen = central_count(params, args%lnm)
        if ( abs(params%sigma_m) < 1e-08_c_double ) then
            args%n_cen = int(p_cen, kind=c_int64_t)
        else
            args%n_cen = binomial_rv(rstate, 1_c_int64_t, p_cen)
        end if
        if ( args%n_cen < 1 ) return ! No central galaxies -> no satellites also 

        ! Satellite galaxy count: this drawn from a poisson distribution with the 
        ! calculated average.
        lam_sat    = satellite_count(params, args%lnm)
        args%n_sat = poisson_rv(rstate, lam_sat)
        
        ! Calculating halo concentration parameter
        args%c = halo_concentration(params, args%lnm, s_spline, ns)

    end subroutine setup_catalog_generation

    subroutine generate_satellites(params, args, rstate, pos, mass) bind(c)
        !! Calculate various parameters for a galaxy catalog generation

        type(hmargs_t), intent(in) :: params 
        !! Halo model parameters

        type(cgargs_t), intent(inout) :: args
        !! Arguments for catalog generation

        integer(c_int64_t), intent(inout) :: rstate
        !! Random number generator state

        real(c_double), intent(out) :: pos(:,:)
        !! Galaxy positions. Must have enough size for holding all the 
        !! galaxies.

        real(c_double), intent(out) :: mass(:)
        !! Galaxy masses.

        integer(c_int64_t) :: i, n_sat
        real(c_double) :: m_halo, r_halo, c_halo, f, r, theta, phi, Ac, k1, k2, p

        if ( args%n_cen < 1 ) return ! No galaxies in this halo

        m_halo = exp(args%lnm) ! Halo mass in Msun
        r_halo = exp(args%lnr) ! Halo radius in Mpc
        c_halo = args%c        ! Halo concentration parameter

        ! Halo has a central galaxy: the position and mass of this galaxy is same
        ! as that of the parent halo.
        pos(1:3,1) = args%pos(1:3) ! in Mpc
        mass(1)    = m_halo        ! in Msun

        if ( args%n_sat < 1 ) return ! No satellite galaxies in this halo

        n_sat = 0 ! This is the actual number of satellite galaxies generated
        do i = 1, args%n_sat
            
            ! Assigning random mass values to the satellite galaxies: 
            ! These masses are drown from a bounded pareto distribution, 
            ! with bounds [m_min, scale_shmf*m_halo] and slope given by 
            ! slope_shmf. 
            ! NOTE: RVs are generated using inverse transform sampling 
            ! (<https://en.wikipedia.org/wiki/Pareto_distribution>)
            p  = -1._c_double / params%slope_shmf
            k1 = exp(params%lnm_min - args%lnm)
            k2 = params%scale_shmf**params%slope_shmf
            if ( params%scale_shmf < k1 ) cycle ! Mass range is non existent: no satellites
            k1 = k1**params%slope_shmf
            f  = ( ( k2 - (k2 - k1) * uniform_rv(rstate) ) / (k1*k2) )**p ! m_sat / m_halo
            
            ! Generating random values corresponding to the distance of 
            ! the galaxy from the halo center. These RVs should follow a 
            ! distribution matching the NFW density profile of the halo. 
            ! Sampling is done using the inverse transformation method. 
            Ac    = uniform_rv(rstate)*( log(1 + c_halo) - c_halo / (1 + c_halo) )
            r     = (r_halo / c_halo) * nfw_c(Ac)
            theta = acos( 2*uniform_rv(rstate) - 1 ) ! -pi to pi
            phi   = 2*pi*uniform_rv(rstate)          !   0 to 2pi
            
            ! Satellite galaxy coordinates x, y, and z in Mpc
            pos(1,i+1) = pos(1,1) + r*sin(theta)*cos(phi)
            pos(2,i+1) = pos(2,1) + r*sin(theta)*sin(phi)
            pos(3,i+1) = pos(3,1) + r*cos(theta)
            
            ! Satellute galaxy mass in Msun
            mass(i+1) = mass(1) * f

            n_sat = n_sat + 1
        end do
        args%n_sat = n_sat ! replace with actual number satellites generated 
        
    end subroutine generate_satellites

    function nfw_c(a) result(res)
        !! Return the value of inverse to the NFW mass function. 

        real(c_double), intent(in) :: a
        real(c_double) :: res, x, p1, p2

        integer, parameter :: f8 = c_double

        ! Using an approximate piecewise rational function for inverting the 
        ! c-A(c) relation. This fit works good upto A(c) ~ 10, with less that  
        ! 10% error. For getting the actual value, one should invert the NFW 
        ! mass function, A(c) = log(1+c) - c/(1+c).  
        x = log10(a)
        if ( a < 1e-03_c_double ) then
            res = 0.50018962_f8*x + 0.15241388_f8
        else if ( a < 2._c_double ) then
            p1  = 2699.40545133_f8 + 2921.07235917_f8*x - 1162.90566455_f8*x**2 
            p2  = 3705.23701996_f8 - 3065.75405505_f8*x -   61.92662277_f8*x**2
            res = p1 / p2
        else
            p1  = 47.25501938_f8 + 32.98237791_f8*x - 38.26172387_f8*x**2  
            p2  = 66.29326139_f8 - 87.37718415_f8*x + 29.86558497_f8*x**2
            res = p1 / p2
        end if
        res = 10._c_double**res
        
    end function nfw_c

! Halo averages:

    function average_halo_density(args, lnma, lnmb, hmf_spline, ns, &
                                  abstol, reltol, maxiter         ) result(res) bind(c)
        !! Return the average halo number density for the given halo mass 
        !! range at current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnma
        !! Lower limit: natural log of halo mass (Msun)
        
        real(c_double), intent(in), value :: lnmb
        !! Upper limit: natural log of halo mass (Msun)

        real(c_double), intent(in) :: hmf_spline(3,ns)
        !! Cubic spline for halo mass function, log(dn/dlnm) as function 
        !! of log(halo mass)

        integer(c_int64_t), intent(in), value :: ns
        !! Size of the mass-function spline data

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double) :: res, err 
        integer(c_int) :: stat
        
        call calc_halo_average(0, lnma, lnmb, args, hmf_spline, ns,    &
                               abstol, reltol, maxiter, res, err, stat &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_halo_density: integral failed to converge'
        
    end function average_halo_density

    function average_galaxy_density(args, lnma, lnmb, hmf_spline, ns, &
                                    abstol, reltol, maxiter         ) result(res) bind(c)
        !! Return the average galaxy number density for the given halo mass 
        !! range at current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnma
        !! Lower limit: natural log of halo mass (Msun)
        
        real(c_double), intent(in), value :: lnmb
        !! Upper limit: natural log of halo mass (Msun)

        real(c_double), intent(in) :: hmf_spline(3,ns)
        !! Cubic spline for halo mass function, log(dn/dlnm) as function 
        !! of log(halo mass)

        integer(c_int64_t), intent(in), value :: ns
        !! Size of the mass-function spline data

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double) :: res, res1, res2, err 
        integer(c_int) :: stat
        
        ! Average central galaxy density
        call calc_halo_average(1, lnma, lnmb, args, hmf_spline, ns,    &
                               abstol, reltol, maxiter, res1, err, stat &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_galaxy_density: integral (c) failed to converge'
        
        ! Average satellite galaxy density
        call calc_halo_average(2, lnma, lnmb, args, hmf_spline, ns,    &
                               abstol, reltol, maxiter, res2, err, stat &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_galaxy_density: integral (s) failed to converge'

        ! Average galaxy density
        res = res1 + res2

    end function average_galaxy_density

    function average_satellite_frac(args, lnma, lnmb, hmf_spline, ns, &
                                    abstol, reltol, maxiter         ) result(res) bind(c)
        !! Return the average satellite fraction for the given halo mass 
        !! range at current redshift.

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in), value :: lnma
        !! Lower limit: natural log of halo mass (Msun)
        
        real(c_double), intent(in), value :: lnmb
        !! Upper limit: natural log of halo mass (Msun)

        real(c_double), intent(in) :: hmf_spline(3,ns)
        !! Cubic spline for halo mass function, log(dn/dlnm) as function 
        !! of log(halo mass)

        integer(c_int64_t), intent(in), value :: ns
        !! Size of the mass-function spline data

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double) :: res, res1, res2, err 
        integer(c_int) :: stat
        
        ! Average central galaxy density
        call calc_halo_average(1, lnma, lnmb, args, hmf_spline, ns,    &
                               abstol, reltol, maxiter, res1, err, stat &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_galaxy_density: integral (c) failed to converge'
        
        ! Average satellite galaxy density
        call calc_halo_average(2, lnma, lnmb, args, hmf_spline, ns,    &
                               abstol, reltol, maxiter, res2, err, stat &
        )
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: average_galaxy_density: integral (s) failed to converge'

        ! Average satellite fraction
        res = res2 / (res1 + res2)
        
    end function average_satellite_frac

    subroutine calc_halo_average(f, a, b, args, hmf_spline, ns, &
                                 abstol, reltol, maxiter, res, err, stat)
        !! Calculate the integral of a scalar function f(x) over the interval [a, b]. 

        integer(c_int), intent(in), value :: f
        !! Function selector (0=no weight, 1=central count, 2=satellite count)

        real(c_double), intent(in), value :: a
        !! Lower limit of integration

        real(c_double), intent(in), value :: b
        !! Upper limit of integration

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in) :: hmf_spline(3,ns)
        !! Mass-function spline

        integer(c_int64_t), intent(in), value :: ns
        !! Size of the mass-function spline

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
        call calc_halo_average2(f, a, b, args, hmf_spline, ns, I0, err0)
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
            call calc_halo_average2(f, xa, xm, args, hmf_spline, ns, I1, err1)
            call int_heap_push(heap, heap_size, heap_capacity, xa, xm, I1, err1) ! Push new interval back
            
            ! Refine on left interval
            call calc_halo_average2(f, xm, xb, args, hmf_spline, ns, I2, err2)
            call int_heap_push(heap, heap_size, heap_capacity, xm, xb, I2, err2) ! Push new interval back
            
            ! Update global sums
            res = res + (I1   + I2   - I0  ) ! replace old interval
            err = err + (err1 + err2 - err0)
            
        end do

        deallocate(heap)    
        
    end subroutine calc_halo_average

    subroutine calc_halo_average2(f, a, b, args, hmf_spline, ns, res, err)
        !! Calculate the integral of a scalar function f(x) over the interval [a, b]. 

        integer(c_int), intent(in), value :: f
        !! Function selector (0=no weight, 1=central count, 2=satellite count)

        real(c_double), intent(in), value :: a
        !! Lower limit of integration

        real(c_double), intent(in), value :: b
        !! Upper limit of integration

        type(hmargs_t), intent(in) :: args
        !! Model parameter values

        real(c_double), intent(in) :: hmf_spline(3,ns)
        !! Mass-function spline

        integer(c_int64_t), intent(in), value :: ns
        !! Size of the mass-function spline

        real(c_double), intent(out) :: res
        !! Value of the integral of f over [a, b]

        real(c_double), intent(out) :: err
        !! Estimate of the error in integration
            
        integer(c_int64_t) :: j
        real(c_double)     :: intg, intk, xval, fval, fval2, scale

        scale = 0.5_c_double * (b - a)
        xval  = a + scale ! log(m)
        fval  = exp( interpolate(xval, ns, hmf_spline) ) ! dn/dlnm
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
            fval = exp( interpolate(xval, ns, hmf_spline) ) ! dn/dlnm
            select case ( f )
            case ( 1 ) ! weight=central galaxy count
                fval = fval * central_count(args, xval)
            case ( 2 ) ! weight=satellite galaxy count
                fval = fval * satellite_count(args, xval)
            ! default=no weight
            end select 

            xval  = a + scale * (1. + K15(1,j)) ! log(m)
            fval2 = exp( interpolate(xval, ns, hmf_spline) ) ! dn/dlnm
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