module hmf_baselib_mod
    !! Values and derived types required for halo mass-function calculations.

    use iso_c_binding
    use variance_mod
    implicit none

    real(c_double), parameter :: delta_sc = 1.6864701998411453_c_double
    !! Overdensity threshold for spherical collapse in EdS universe

    real(c_double), parameter :: critical_density_const = 2.775e+11_c_double
    !! Constant part of the present critical density

    real(c_double), parameter, private :: pi  = 3.141592653589793_c_double
    !! Pi

    type, bind(c) :: hmfargs_t
        !! A struct containing values of various arguments for halo 
        !! mass-function calculation routines.

        real(c_double) :: z
        !! Redshift

        real(c_double) :: H0
        !! Hubble parameter
        
        real(c_double) :: Om0
        !! Total matter density parameter

        real(c_double) :: Delta_m
        !! Matter overdensity w.r.to mean background density

        real(c_double) :: lnm
        !! Natural log of halo mass in Msun

        real(c_double) :: s
        !! Matter variance corresponding to halo mass

        real(c_double) :: dlnsdlnm 
        !! LOg derivative of matter variance w.r.to halo mass

        real(c_double) :: rho_m
        !! Total matter density at redshift 0 (unit: Msun/Mpc^3)

        real(c_double) :: param(6)
        !! Model parameters
        
    end type
    
contains

    subroutine setup_hmf_calculation(args, filt, pktab, size, cls)
        !! Calculate related quantities for halo mass-function calculation.
        !!
        !! NOTE: setup method for the specific mass-function model (if any) should be 
        !! run after this.

        type(hmfargs_t), intent(inout) :: args

        integer(c_int), intent(in), value :: filt
        !! Code of filter function (0=tophat, 1=gaussian)

        real(c_double), intent(in) :: pktab(cls, size)
        !! Precalculated power spectrum table. The columns should be
        !! 1=Nodes of integration (natural log of k in 1/Mpc), 
        !! 2=Value of natural log of power spectrum, 
        !! 3=Weights for integration.
        
        integer(c_int64_t), intent(in), value :: size
        !! Size of the power spectrum table: must be a multiple of 7.

        integer(c_int), intent(in), value :: cls
        !! Columns of the power spectrum table: must be 2 or 3. If 3, use the 
        !! the last column as weights. Otherwise, Simpson's rule is used. For 
        !! other values, return value will be NaN.

        real(c_double) :: lnr, rho_h, s2

        args%rho_m = ( critical_density_const * args%H0**2 ) ! matter density at z=0 in Msun/Mpc^3 

        ! Lagrangian radius (r) corresponding to halo mass
        rho_h = args%rho_m ! * args%Delta_m
        lnr   = ( args%lnm + log(3._c_double / (4*pi) / rho_h ) ) / 3._c_double ! r in Mpc

        ! Matter variance
        s2     = variance(lnr, 0, 0, filt, pktab, size, cls)
        args%s = sqrt(s2)
        
        ! Log derivative of matter variance
        args%dlnsdlnm = variance(lnr, 0, 1, filt, pktab, size, cls) ! ds2/dr
        args%dlnsdlnm = args%dlnsdlnm * ( exp(lnr) / s2 / 6._c_double )
        
    end subroutine setup_hmf_calculation

    function convert_fs_to_hmf(args, fs, target_code) result(res)
        !! Calculate the value of halo mass-function given f(s) and other values.

        type(hmfargs_t), intent(in) :: args
        !! Arguments

        real(c_double), intent(in), value :: fs
        !! Value of mass-function f(s)

        integer(c_int), intent(in), value :: target_code
        !! Code for output value (0=dn/dm, 1=dn/dlnm, 2=dn/dlog10m)
        
        real(c_double) :: res, m

        m = exp(args%lnm) ! Halo mass in Msun
        
        res = fs * abs(args%dlnsdlnm) * args%rho_m / m ! dn/dlnm in Mpc^-3
        select case ( target_code )
        case ( 2 )
            res = res / log(10._c_double) ! dn/dlog10m
        case ( 0 ) 
            res = res / m ! dn/dm
        end select

    end function convert_fs_to_hmf
    
end module hmf_baselib_mod