module growthfactor
    !!
    !! Calculation of linear growth factor in a w0-wa CDM cosmology model.
    !!
    use iso_c_binding
    use quadutils
    implicit none

    private :: integrand

    type, public, bind(c) :: lgargs_t
    !! A struct containing values of various arguments for growth factor 
    !! calculation routines.
        
        real(c_double) :: Om0 
        !! Total matter density parameter

        real(c_double) :: Ode0
        !! Dark energy density parameter

        real(c_double) :: w0
        !! Present value of the dark energy equation of state, w
        
        real(c_double) :: wa
        !! A measure of how w evolves with time

        real(c_double) :: abstol
        !! Absolute tolerance for checking convergence in integration

        real(c_double) :: reltol
        !! Relative tolerance for checking convergence in integration

        integer(c_int64_t) :: maxiter
        !! Maximum number of iterations for calculating integral

    end type
    
contains

    function integrand(a, argsp) result(res)
        !! Integrand for growth factor calculations

        real(c_double), value :: a
        !! Scale factor 

        type(c_ptr), value :: argsp
        !! Model parameters (w0-wa CDM model)
        
        real(c_double) :: res
        !! Value of the function

        type(lgargs_t), pointer :: cm ! w0-wa CDM model
        real(c_double) :: p, q ! Related to drak energy density

        ! Get pointer to args 
        call c_f_pointer(argsp, cm)

        ! Calculating Hubble function, E^2(a)
        p   = 3*cm%wa
        q   = 3*(1 + cm%w0 + cm%wa) 
        res = ( cm%Om0 / a + (1 - cm%Om0 - cm%Ode0) ) / a**2 ! Matter + Curvature
        res = res + cm%Ode0 * exp( p*(a - 1) ) * a**(-q)     ! w0-wa dark energy
        res = sqrt(res)

        ! Integrand
        res = (res * a)**(-3)

    end function integrand

    function linear_growth(z, nu, args) result(res) bind(c)
        !! Calculate the value of linear growth factor in a w0-wa CDM model
        !! cosmology.

        real(c_double), intent(in), value :: z
        !! Redshift

        integer(c_int), intent(in), value :: nu
        !! Return value code (0=growth factor, 1=growth rate)

        type(lgargs_t), intent(in), target :: args
        !! Values for model and control parameters
        
        real(c_double) :: res
        !! Return value

        real(c_double), parameter :: a_start = 0.0_c_double

        real(c_double) :: abstol, reltol, a, ym, yk, yde, p, q, y, err
        integer(c_int) :: stat
        integer(c_int64_t) :: maxiter

        ! Set values for optional control parameters
        abstol  = max(args%abstol , 1e-08_c_double)
        reltol  = max(args%reltol , 1e-08_c_double)
        maxiter = max(args%maxiter, 50_c_int64_t  )
        
        a = 1.0_c_double / (z + 1.0_c_double) ! scale factor
        call integrate(integrand, a_start, a, c_loc(args), &
                       abstol, reltol, maxiter,            &
                       res, err, stat                      &
        )

                       
        ! Calculating Hubble function, E^2(a)
        ym  = args%Om0 / a**3                     ! matter
        yk  = ( 1 - args%Om0 - args%Ode0 ) / a**2 ! curvature
        p   = 3*args%wa * (a - 1)
        q   = 3*(1 + args%w0 + args%wa)  
        yde = args%Ode0 * exp(p) * a**(-q)        ! w0-wa dark energy
        y   = ym + yk + yde

        if ( nu /= 0 ) then
            ! Linear growth rate: first log derivative
            res = 1.0_c_double / ( (a*y)**3 * res )
            p   = p + q - 3*args%wa 
            y   = -( 3*ym + 2*yk + p*yde ) / (2*y) ! log derivative of E(a)
            res = res + y
        else
            ! Linear growth factor
            res = 2.5_c_double*args%Om0 * sqrt(y) * res
        end if
        
    end function linear_growth
    
end module growthfactor