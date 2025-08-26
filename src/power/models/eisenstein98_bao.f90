module eisenstein98_bao_mod
    !! Eisenstein & Hu (1998) model for linear matter power spectrum inclding 
    !! baryon acoustic oscillations.

    use iso_c_binding
    use pm_ingredients_mod
    implicit none

    private
    public :: init_eisenstein98_bao, ps_eisenstein98_bao

contains

    subroutine init_eisenstein98_bao(args) bind(c)
        !! Initialize parameters related to Eisenstein & Hu matter power spectrum
        !! including baryon acoustic oscillations (BAO).

        type(psargs_t), intent(inout) :: args

        real(c_double) :: Omh2, Obh2, fb, fc, c1, c2, k_eq, R_eq, R_d, y, Gy

        Omh2 = args%Omh2
        Obh2 = args%Obh2
        fb   = Obh2 / Omh2                  ! Baryon fraction
        fc   = 1 - fb - (args%Onuh2 / Omh2) ! CDM fraction

        ! Redshift at matter-radiation equality (eqn. 1)
        args%z_eq = 2.5e+04_dp * Omh2 / args%theta**4

        ! Redshift at drag epoch (eqn. 2)
        c1  = 0.313_dp*(1 + 0.607_dp*Omh2**0.674_dp) / Omh2**0.419_dp
        c2  = 0.238_dp*Omh2**0.223_dp
        args%z_d = 1291.0_dp*(Omh2**0.251_dp)*(1 + c1*Obh2**c2) / (1 + 0.659_dp*Omh2**0.828_dp)

        ! Scale of particle horizon at z_eq
        k_eq = 7.46e-02_dp*Omh2 / args%theta**2

        ! Ratio of baryon - photon momentum density (eqn. 5)
        R_eq = 31.5_dp*Obh2*args%theta**(-4) * (args%z_eq / 1.0e+03_dp)**(-1) ! at z_eq
        R_d  = 31.5_dp*Obh2*args%theta**(-4) * (args%z_d  / 1.0e+03_dp)**(-1) ! at z_d

        ! Sound horizon (eqn. 26)
        args%s = (2._dp / 3._dp / k_eq)    &
                    * sqrt(6._dp / R_eq)   &
                    * log((sqrt(1 + R_d) + sqrt(R_d + R_eq)) / (1 + sqrt(R_eq)))

        ! Silk damping scale (eqn. 7)
        args%k_silk = 1.6_dp*Obh2**0.52_dp * Omh2**0.73_dp * (1 + (10.4_dp*Omh2)**(-0.95_dp))

        ! Parameter alpha_c, Eqn. 11
        c1 = (46.9_dp*Omh2)**0.670_dp * (1._dp + (32.1_dp*Omh2)**(-0.532_dp))
        c2 = (12.0_dp*Omh2)**0.424_dp * (1._dp + (45.0_dp*Omh2)**(-0.582_dp))
        args%param(1) = c1**(-fb) * c2**(-fb**3)

        ! Parameter beta_c, Eqn. 12
        c1 = 0.944_dp*(1._dp + (458._dp*Omh2)**(-0.708_dp))**(-1)
        c2 = (0.395_dp*Omh2)**(-0.0266_dp)
        args%param(2) = (1._dp + c1*(fc**c2 - 1._dp))**(-1)

        ! Parameter alpha_b, Eqn. 14-15
        y  = (1 + args%z_eq) / (1 + args%z_d)
        Gy = y*( -6*sqrt(1 + y) + (2 + 3*y) * log( (sqrt(1 + y) + 1) / (sqrt(1 + y) - 1) ) )
        args%param(3) = 2.07_dp*k_eq * args%s * (1 + R_d)**(-0.75_dp) * Gy

        ! Parameter beta_b, Eqn. 24
        args%param(4) = 0.5_dp + fb + (3 - 2*fb)*sqrt((17.2_dp*Omh2)**2 + 1)

        ! Parameter beta_node, Eqn. 23
        args%param(5) = 8.41_dp*Omh2**0.435_dp

    end subroutine init_eisenstein98_bao

    function ps_eisenstein98_bao(lnk, args) result(retval) bind(c)
        !! Return the value of Eisenstein & Hu matter power spectrum including baryon
        !! acoustic oscillations (BAO).
        
        real(c_double), intent(in), value :: lnk
        !! Natural log of wavenumber in 1/Mpc

        type(psargs_t), intent(in) :: args
        !! Other arguments

        real(c_double) :: retval
        !! Natural log of matter power spectrum value

        real(c_double) :: alpha_c, alpha_b, beta_c, beta_b, beta_node, dplus_z,   &
                        fc, fb, s, ks, q, k, t0b, t0ab, t2, t3, t4, f, st

        alpha_c   = args%param(1)
        alpha_b   = args%param(2)
        beta_c    = args%param(3) 
        beta_b    = args%param(4)
        beta_node = args%param(5)
        s         = args%s
        fb        = args%Obh2 / args%Omh2             ! Baryon fraction
        fc        = 1 - fb - (args%Onuh2 / args%Omh2) ! CDM fraction
        
        ! Linear growth factor, normalized so that value = 1 at redshift 0
        dplus_z = args%dplus_z / args%dplus_0 

        k  = exp(lnk) ! Mpc^-1
        q  = k * ( args%theta**2 / args%Omh2 ) 
        ks = k*s 

        ! Transfer function: CDM part (Eqn. 17-20)
        f    = 1._dp / (1 + (ks / 5.4_dp)**4)  ! Eqn. 18
        t2   = log( e + 1.8_dp*beta_c*q )
        t3   = 14.2_dp + 386.0_dp / (1 + 69.9_dp*q**1.08_dp)            ! Eqn. 20 with alpha_c=1
        t4   = 14.2_dp / alpha_c + 386.0_dp / (1 + 69.9_dp*q**1.08_dp ) ! Eqn. 20 with alpha_c
        t0b  = t2 / (t2 + t3*q**2)           ! Eqn. 19 
        t0ab = t2 / (t2 + t4*q**2)           ! Eqn. 19 
        retval = fc*( f*t0b + (1 - f)*t0ab ) ! Eqn. 17

        ! Transfer function: baryon part (Eqn. 21)
        st   = s / (1 + (beta_node/ks)**3)**(1._dp/3._dp)
        f    = sin(k*st) / (k*st) ! Spherical Bessel function, j0
        t2   = log( e + 1.8_dp*q )
        t0ab = t2 / (t2 + t3*q**2)
        t0b  = t0ab / (1 + (ks/5.2_dp)**2)     &
                + alpha_b / (1 + (beta_b/ks)**3) * exp(-(k/args%k_silk)**1.4_dp)
        retval = retval  + fb * t0b * f

        ! Interpolation using growth factor
        retval = dplus_z * retval 

        ! Linear matter power spectrum
        retval = (args%sigma8**2 * args%norm) * retval**2 * k**( args%ns )
        retval = log(retval)
        
    end function ps_eisenstein98_bao
    
end module eisenstein98_bao_mod