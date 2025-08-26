module eisenstein98_mnu_mod
    !! Eisenstein & Hu (1998) model for linear matter power spectrum
    !! including massive neutrinos.
    
    use iso_c_binding
    use pm_ingredients_mod
    implicit none

    private
    public :: init_eisenstein98_mnu, ps_eisenstein98_mnu
    
contains

    subroutine init_eisenstein98_mnu(args) bind(c)
        !! Initialize parameters related to Eisenstein & Hu matter power spectrum
        !! including massive neutrino.

        type(psargs_t), intent(inout) :: args

        real(c_double) :: Omh2, Obh2, Onuh2, Nnu, fb, fnu, fnb, fc, fcb, pc, pcb, &
                        c1, c2, yd

        Omh2  = args%Omh2
        Obh2  = args%Obh2
        Onuh2 = args%Onuh2
        Nnu   = args%Nnu
        fb    = Obh2  / Omh2 ! Baryon fraction     
        fnu   = Onuh2 / Omh2 ! Massive neutrino fraction
        fnb   = fnu + fb     ! Massive neutrino + baryon    
        fc    = 1._dp - fnb  ! CDM fraction    
        fcb   = fc + fb      ! CDM + baryon fraction

        ! Redshift at matter-radiation equality (eqn. 1)
        args%z_eq = 2.5e+04_dp * Omh2 / args%theta**4 

        ! Redshift at drag epoch (eqn. 2)
        c1  = 0.313_dp*(1 + 0.607_dp*Omh2**0.674_dp) / Omh2**0.419_dp
        c2  = 0.238_dp*Omh2**0.223_dp
        args%z_d = 1291.0_dp*(Omh2**0.251_dp)*(1 + c1*Obh2**c2) / (1 + 0.659_dp*Omh2**0.828_dp)

        ! Sound horizon (eqn. 26)
        args%s = 44.5_dp*log( 9.83_dp/Omh2 ) / sqrt( 1._dp + 10._dp*Obh2**0.75_dp ) 

        ! Eqn. 14 in EH98 paper
        pc  = 0.25_dp*( 5._dp - sqrt( 1._dp + 24.0_dp*fc  ) ) 
        pcb = 0.25_dp*( 5._dp - sqrt( 1._dp + 24.0_dp*fcb ) )

        ! Small-scale suppression, alpha_nu (Eqn. 15)
        yd  = (1 + args%z_eq) / (1 + args%z_d) ! eqn. 3
        args%param(1) = (fc / fcb) * (5 - 2*(pc + pcb)) / (5 - 4*pcb)   &                                                 
            * (1 - 0.533_dp*fnb + 0.126_dp*fnb**3)                      &
            / (1 - 0.193_dp*sqrt(fnu*Nnu) + 0.169_dp*fnu*Nnu**0.2_dp)   &
            * (1 + yd)**(pcb - pc)                                      &                                        
            * (1 + 0.5_dp*(pc - pcb) * (1 + 1._dp / (3 - 4*pc) / (7 - 4*pcb)) * (1 + yd)**(-1))
        
        ! Eqn. 21, beta_c
        args%param(2) = (1 - 0.949_dp*fnb)**(-1)

        ! Parameter related to growth factor suppression (y_fs / q^2) (eqn. 14)
        args%param(3) = 17.2_dp*fnu*( 1 + 0.488_dp*fnu**(-7._dp/6._dp) ) * (Nnu/fnu)**2

        ! Constant part of B_k in Eqn. 22
        args%param(4) = 1.2_dp*fnu**0.64_dp * Nnu**(0.3_dp + 0.6_dp*fnu)

    end subroutine init_eisenstein98_mnu

    function ps_eisenstein98_mnu(lnk, args) result(retval) bind(c)
        !! Return the value of Eisenstein & Hu matter power spectrum including massive
        !! neutrinos.
        
        real(c_double), intent(in), value :: lnk
        !! Natural log of wavenumber in 1/Mpc

        type(psargs_t), intent(in) :: args
        !! Other arguments

        real(c_double) :: retval
        !! Natural log of matter power spectrum value

        real(c_double) :: Omh2, anu, dplus_z, dplus_c, beta_c, fnu, fcb, pcb,   &
                        Bconst, yfs_q2, s, q, k, t2, t3

        Omh2   = args%Omh2
        anu    = sqrt(args%param(1))
        beta_c = args%param(2)
        yfs_q2 = args%param(3)
        Bconst = args%param(4)
        s      = args%s
        fnu    = args%Onuh2 / Omh2 ! Massive neutrino fraction
        fcb    = 1 - fnu           ! CDM + baryon fraction
        pcb    = 0.25_dp*( 5._dp - sqrt( 1._dp + 24.0_dp*fcb ) )
        
        k = exp(lnk) ! wavenumber in Mpc^-1
        q = k * ( args%theta**2 / Omh2 ) ! dimension-less wavenumber

        ! Linear growth factor
        dplus_c = 2.5_dp*( Omh2 / args%h**2 )*(args%z_eq + 1) ! Used as normalization factor for growth
        dplus_z = dplus_c * args%dplus_z

        ! Suppressed growth factor
        t2  = yfs_q2 * q**2 ! y_fs in eqn. 14
        if ( args%include_nu /= 0 ) then 
            ! Dcbnu (eqn., 12), including neutrino
            t2 = ( fcb**(0.7_dp/pcb) + dplus_z / (1 + t2) )**(pcb/0.7_dp) * dplus_z**(1 - pcb) 
        else 
            ! Dcb (eqn., 13), not incl. neutrino
            t2 = ( (1 + dplus_z) / (1 + t2) )**(pcb/0.7_dp) * dplus_z**(1 - pcb) 
        end if
        retval = t2 / dplus_z

        ! Master function T_master (Eqn. 22-24)
        q  = 3.92_dp*q*sqrt(args%Nnu) / fnu                 ! q_nu in Eqn. 23
        t2 = 1 + Bconst / (q**(-1.6_dp) + q**0.8_dp) ! B_k in Eqn. 22
        retval = retval * t2   

        ! Transfer function T_sup (Eqn. 17-20)
        t2 = Omh2*( anu + (1 - anu)/(1 + (0.43*k*s)**4) ) ! Shape parameter, Gamma_eff (Eqn. 16)
        q  = k * args%theta**2 / t2                       ! Effective wavenumber (Eqn. 17)
        t2 = log( e + 1.84_dp*beta_c * anu * q )          ! L (Eqn. 19)
        t3 = 14.4_dp + 325._dp / (1 + 60.5_dp*q**1.11_dp) ! C (Eqn. 20)
        retval = t2 / (t2 + t3*q**2)                      ! Eqn. 18

        ! Linear interpolation using growth factor
        dplus_z = args%dplus_z / args%dplus_0 ! normalized so that value = 1 at redshift 0
        retval  = retval * dplus_z  

        ! Linear matter power spectrum
        retval = (args%sigma8**2 * args%norm) * retval**2 * k**( args%ns )
        retval = log(retval)
        
    end function ps_eisenstein98_mnu
    
end module eisenstein98_mnu_mod