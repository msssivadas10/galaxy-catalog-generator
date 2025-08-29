module zfunctions_mod
    !! Calculation of some functions of redshiftin a w0-wa CDM cosmology model, 
    !! such as the linear growth factor.

    use iso_fortran_env, only: stderr => error_unit
    use iso_c_binding
    use integrate_mod
    implicit none

    private 
    public :: linear_growth

    type, public, bind(c) :: zfargs_t
        !! A struct containing values of various arguments for the redshift 
        !! function calculation routines.
        
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

! Linear Growth Factor:

    function growth_integrand(a, args) result(res)
        !! Integrand for growth factor calculations: `(a*E(a))**-1.5`

        real(c_double), value :: a
        !! Scale factor 

        type(zfargs_t) :: args
        !! Model parameters (w0-wa CDM model)
        
        real(c_double) :: res
        !! Value of the function

        real(c_double) :: p, q, Ok0

        Ok0 = (1 - args%Om0 - args%Ode0)
        p   = 3*args%wa
        q   = 3*(1 + args%w0 + args%wa) 
        res = a / ( args%Om0 + Ok0 * a + args%Ode0 * exp( p*(a-1) ) * a**(3-q) )
        res = res**(1.5_c_double)

    end function growth_integrand

    function linear_growth(z, nu, args) result(res) bind(c)
        !! Calculate the value of linear growth factor in a w0-wa CDM model
        !! cosmology.

        real(c_double), intent(in), value :: z
        !! Redshift

        integer(c_int), intent(in), value :: nu
        !! Return value code (0=growth factor, 1=growth rate)

        type(zfargs_t), intent(in), target :: args
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
        call integrate(growth_integrand, a_start, a, args, &
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

! Helper functions:

    subroutine integrate(f, a, b, args, abstol, reltol, maxiter, res, err, stat)
        !! Calculate the integral of a scalar function f(x) over the interval [a, b]. 

        interface
            function f(x, args_) result(y) bind(c)
                !! Function to integrate
                import :: c_double, zfargs_t
                real(c_double), value :: x
                type(zfargs_t) :: args_ 
                real(c_double) :: y
            end function
        end interface

        real(c_double), intent(in), value :: a
        !! Lower limit of integration

        real(c_double), intent(in), value :: b
        !! Upper limit of integration

        type(zfargs_t) :: args
        !! Other arguments to pass to the function

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

        integer(c_int64_t) :: iter, j
        real(c_double)     :: xval(15), fval(15)
        real(c_double)     :: xa, xb, xm, I0, I1, I2, err0, err1, err2
        integer(c_int64_t) :: heap_size, heap_capacity
        real(c_double), allocatable :: heap(:, :)

        heap_size     = 0
        heap_capacity = 10*maxiter ! Heap capacity
        allocate( heap(4, heap_capacity) )

        ! Initial evaluation
        call G7K15_nodes(a, b, xval)             ! Get integration nodes 
        fval = [( f(xval(j), args), j = 1, 15 )] ! Evaluate f(x) at nodes
        call G7K15(fval, a, b, I0, err0)
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
            call G7K15_nodes(xa, xm, xval)           ! Get integration nodes 
            fval = [( f(xval(j), args), j = 1, 15 )] ! Evaluate f(x) at nodes
            call G7K15(fval, xa, xm, I1, err1)
            call int_heap_push(heap, heap_size, heap_capacity, xa, xm, I1, err1) ! Push new interval back
            
            ! Refine on left interval
            call G7K15_nodes(xm, xb, xval)           ! Get integration nodes 
            fval = [( f(xval(j), args), j = 1, 15 )] ! Evaluate f(x) at nodes
            call G7K15(fval, xm, xb, I2, err2)
            call int_heap_push(heap, heap_size, heap_capacity, xm, xb, I2, err2) ! Push new interval back
            
            ! Update global sums
            res = res + (I1   + I2   - I0  ) ! replace old interval
            err = err + (err1 + err2 - err0)
            
        end do
        
        deallocate(heap)
        
    end subroutine integrate
    
end module zfunctions_mod