module pm_ingredients_mod
    !! Values and derived types required for power spectrum calculations.

    use iso_c_binding
    implicit none

    integer, parameter :: dp = c_double

    real(c_double), parameter :: e = 2.718281828459045_dp 
    !! Base of natural logarithm, e
    
    type, bind(c) :: psargs_t
        !! A struct containing values of various arguments and parameters required
        !! for matter power spectrum calculations.

        ! -- Parameters set by the user -- 

        real(c_double) :: z
        !! Redshift at which the power spectrum is evaluated 

        real(c_double) :: h 
        !! Hubble parameter (in unit of 100 km/sec/Mpc)

        real(c_double) :: Omh2
        !! Matter density parameter 

        real(c_double) :: Obh2
        !! Baryon matter density parameter 

        real(c_double) :: Onuh2
        !! Massive neutrino density parameter 

        real(c_double) :: Nnu
        !! Number of massive neutrino species
        
        real(c_double) :: ns
        !! Index of the initial power spectrum

        real(c_double) :: sigma8
        !! Matter variance at scale 8 Mpc/h: to normalize power spectrum

        real(c_double) :: theta
        !! Temperature of the cosmic microwave background (in unit of 2.7 K)

        real(c_double) :: dplus_z
        !! Linear growth at this redshift 
        
        real(c_double) :: dplus_0
        !! Linear growth at redshift 0 (used for normalising the growth factor)
        
        real(c_double) :: norm
        !! Power spectrum normalization for sigma8=1
        
        integer(c_int) :: include_nu
        !! Flag indicating if to include neutrino in the transfer function

        ! -- Parameters to be calculated by init subroutines -- 

        real(c_double) :: z_eq
        !! Redshift corresponding to the matter - radiation equality epoch

        real(c_double) :: z_d
        !! Redshift corresponding to the drag epoch
        
        real(c_double) :: s
        !! Sound horizon in Mpc (see Eqn 26 in the EH paper)

        real(c_double) :: k_silk
        !! Silk damping scale  (see Eqn 7 in the EH paper)

        real(c_double) :: param(5)
        !! Various parameters for the model 
        
    end type

end module pm_ingredients_mod