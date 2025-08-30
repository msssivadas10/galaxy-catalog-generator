module halo_model_mod
    !! Halo occupation distributions and related calculations.

    use iso_c_binding
    implicit none

    ! private

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

    end type
    
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
        !! its mass. This will be a power law specified by the parameters 
        !! `(m0, m1, alpha)`.

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
    
end module halo_model_mod