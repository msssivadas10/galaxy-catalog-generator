module growthfactor_mod
    !! Calculation of linear growth factor in a w0-wa CDM cosmology model.

    use iso_fortran_env, only: stderr => error_unit
    use iso_c_binding
    use integrate_mod
    implicit none

    private 
    public :: linear_growth, solve_growth_ode

    type, public, bind(c) :: lgargs_t
        !! A struct containing values of various arguments for growth 
        !! factor calculation routines.
        
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
        !! Integrand for growth factor calculations: `(a*E(a))**-1.5`

        real(c_double), value :: a
        !! Scale factor 

        type(c_ptr), value :: argsp
        !! Model parameters (w0-wa CDM model)
        
        real(c_double) :: res
        !! Value of the function

        type(lgargs_t), pointer :: cm ! w0-wa CDM model
        real(c_double) :: p, q, Ok0

        ! Get pointer to args 
        call c_f_pointer(argsp, cm)

        Ok0 = (1 - cm%Om0 - cm%Ode0)
        p   = 3*cm%wa
        q   = 3*(1 + cm%w0 + cm%wa) 
        res = a / ( cm%Om0 + Ok0 * a + cm%Ode0 * exp( p*(a-1) ) * a**(3-q) )
        res = res**(1.5_c_double)

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

        real(c_double), parameter :: a_start = 1.e-08_c_double

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
        if ( stat /= 0 ) &
            write(stderr,'(a)') 'warning: linear_growth: integral failed to converge'
                       
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

    subroutine growth_ode(args, a, y)
        
        type(lgargs_t), intent(in) :: args
        !! Values for model and control parameters

        real(c_double), intent(in) :: a
        !! Scale factor
        
        real(c_double), intent(inout) :: y(2)
        !! Current values (replaced by derivative values)

        real(c_double) :: p, q, ym, yk, yde, y1, y2, y3, ydot3

        ! Hubble function E(a) and its derivative
        ym    = args%Om0 / a**3                     ! matter 
        yk    = ( 1 - args%Om0 - args%Ode0 ) / a**2 ! curvature
        p     = 3*args%wa * (a - 1)
        q     = 3*(1 + args%w0 + args%wa)  
        yde   = args%Ode0 * exp(p) * a**(-q)        ! w0-wa dark energy
        y3    = ym + yk + yde                       ! E(a)^2
        p     = p + q - 3*args%wa 
        ydot3 = -( 3*ym + 2*yk + p*yde ) / (2*a*y3) ! dlnE(a)/da

        ! Growth factor ODE
        y1   = y(1)
        y2   = y(2)
        y(1) = y2
        y(2) = -(3.d0 / a + ydot3) + (1.5d0*ym / a**2) * (y1 / y3)
        
    end subroutine growth_ode

    subroutine solve_growth_ode(args, steps, z, dz) bind(c)
        !! Calculate the value of linear growth factor in a w0-wa CDM model
        !! cosmology. This solves the ODE for growth factor using RK4 method
        !! and returns the values in scale factor range [0, 1].
        
        type(lgargs_t), intent(in) :: args
        !! Values for model and control parameters

        integer(c_int64_t), intent(in), value :: steps
        !! Number of steps

        real(c_double), intent(out) :: z(steps)
        !! Redshift values

        real(c_double), intent(out) :: dz(2, steps)
        !! Values of linear growth factor and linear growth rate

        real(c_double), dimension(2) :: y, k1, k2, k3, k4
        real(c_double)     :: a, da, da_2, da_6
        integer(c_int64_t) :: step

        ! Calculate stepsize
        da   = 1._c_double / steps
        da_2 = da / 2._c_double
        da_6 = da / 6._c_double 

        ! Initial values: using D(a) ~ a for smaller a
        y = [ 0._c_double, 1._c_double ]
        
        a = da 
        do step = 1, steps
            ! ODE solution using RK4
            
            k1(:) = y(:) 
            call growth_ode(args, a, k1)

            a  = a + da_2
            k2 = y + k1 * da_2
            call growth_ode(args, a, k2)
            
            k3 = y + k2 * da_2
            call growth_ode(args, a, k3)
            
            a  = a + da_2
            k4 = y + k3 * da
            call growth_ode(args, a, k4)

            dz(:, step) = y + (k1 + 2*k2 + 2*k3 + k4) * da_6
            z(step)     = 1._c_double / a - 1

        end do
        
    end subroutine solve_growth_ode
    
end module growthfactor_mod