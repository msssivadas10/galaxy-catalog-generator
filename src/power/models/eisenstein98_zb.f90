module eisenstein98_zb_mod
    !! Eisenstein & Hu (1998) model for linear matter power spectrum for 
    !! zero baryon case.

    use iso_c_binding
    use pm_ingredients_mod
    implicit none

    private
    public :: init_eisenstein98_zb, ps_eisenstein98_zb
    
contains

    subroutine init_eisenstein98_zb(args) bind(c)
        !! Initialize parameters related to Eisenstein & Hu matter power spectrum
        !! for zero baryon.

        type(psargs_t), intent(inout) :: args

        real(c_double) :: Omh2, Obh2, theta, c1, c2, fb

        Omh2  = args%Omh2
        Obh2  = args%Obh2
        theta = args%theta

        ! Redshift at matter-radiation equality (Eqn. 1)
        args%z_eq = 2.5e+04_dp * Omh2 / theta**4

        ! Redshift at drag epoch (Eqn. 2)
        c1  = 0.313_dp*(1 + 0.607_dp*Omh2**0.674_dp) / Omh2**0.419_dp
        c2  = 0.238_dp*Omh2**0.223_dp
        args%z_d = 1291.0_dp*(Omh2**0.251_dp)*(1 + c1*Obh2**c2) / (1 + 0.659_dp*Omh2**0.828_dp)

        ! Sound horizon (Eqn. 26)
        args%s = 44.5_dp*log( 9.83_dp/Omh2 ) / sqrt( 1 + 10._dp*Obh2**0.75_dp )

        ! Parameter alpha_Gamma, Eqn. 31
        fb = Obh2 / Omh2
        args%param(1) = 1._dp - 0.328_dp*log( 431*Omh2 ) * fb + 0.38_dp*log( 22.3_dp*Omh2 ) * fb**2

    end subroutine init_eisenstein98_zb

    function ps_eisenstein98_zb(lnk, args) result(retval) bind(c)
        !! Return the value of Eisenstein & Hu matter power spectrum for zero baryon.
        
        real(c_double), intent(in), value :: lnk
        !! Natural log of wavenumber in 1/Mpc

        type(psargs_t), intent(in) :: args
        !! Other arguments

        real(c_double) :: retval 
        !! Natural log of matter power spectrum value

        real(c_double) :: Omh2, theta, alpha_g, gamma_eff, dplus_z, s, q, k, t2, t3

        Omh2    = args%Omh2
        theta   = args%theta 
        alpha_g = args%param(1)
        s       = args%s
        
        ! Linear growth factor, normalized so that value = 1 at redshift 0
        dplus_z = args%dplus_z / args%dplus_0 

        k = exp(lnk) ! Mpc^-1

        ! Shape parameter, Eqn 30
        gamma_eff = Omh2 * ( alpha_g + ( 1 - alpha_g ) / ( 1 + ( 0.43_dp*k*s )**4 ) ) 
        q = k * ( theta**2 / gamma_eff ) 
        
        ! Transfer function
        t2 = log( 2*e + 1.8_dp*q )
        t3 = 14.2_dp + 731.0_dp / ( 1 + 62.5_dp*q )
        retval = t2 / (t2 + t3*q**2)

        ! Interpolation using growth factor
        retval = dplus_z * retval 

        ! Linear matter power spectrum
        retval = (args%sigma8**2 * args%norm) * retval**2 * k**( args%ns )
        retval = log(retval)
        
    end function ps_eisenstein98_zb
    
end module eisenstein98_zb_mod