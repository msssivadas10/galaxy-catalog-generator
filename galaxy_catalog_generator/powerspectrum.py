
import numpy as np
from scipy.integrate import quad
from functools import reduce
from typing import TypeVar, Literal
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from astropy.cosmology import FLRW

_T = TypeVar('_T')

def linearGrowth(
        z       : _T, 
        cosmo   : FLRW, 
        nu      : Literal[0, 1] = 0, 
        stepsize: float         = 0.01, 
    ) -> _T:
    r"""
    Calculate the linear growth factor, or its logarithmic derivative, according to the specified 
    cosmology model.

    Parameters
    ----------
    z : array_like
        Redshift - values must be greater than -1.
    
    cosmo : astropy.cosmology.FLRW
        Cosmology model. 
        
        .. note:: 
        
            Can be any object which provides the hubble parameter evolution as function of redshift,  
            in the name ``efunc``, with signature ``(float) -> float``.

    nu : {0, 1}, optional
        Return the linear growth rate (first log-derivative of the growth factor) for non-zero value.

    stepsize : float, default=0.01
        Stepsize used to numerically calculate the derivative using finite difference.

    Returns
    -------
    retval : array_like
        Linear growth factor or its first log-derivative.

    Usage
    -----

    >>> from astropy.cosmology import FlatLambdaCDM
    >>> cm = FlatLambdaCDM(H0 = 70., Om0 = 0.3, Ob0 = 0.05, Tcmb0 = 2.725)
    >>> z  = [0., 1.]
    >>> linearGrowth(z, cm) # linear growth factor
    array([1.03798112, 0.63471574])
    >>> linearGrowth(z, cm, nu = 1) # linear growth rate
    array([0.51323816, 0.87038474])

    """
    
    a  = 1. / ( np.asarray(z, dtype=np.float64) + 1. ) # scale factor
    Ez = cosmo.efunc(z) # hubble parameter function

    # -- Linear growth factor calculations:

    retval = np.zeros_like(a)
    
    # For scale factor values close to 0, growth factor is approximately same as the scale factor. 
    # This is used to skip integration for smaller scale factor values...
    mask  = ( np.abs(a) > 1e-06 )
    retval[~mask] = a[~mask] 

    # Calculate the integral expression for other values of scale factor...
    quadVectorized = np.vectorize(quad, otypes = [np.float64, np.float64]) # vectorized version of quad
    integ, _abserr = quadVectorized(
        lambda a: ( a * cosmo.efunc( 1./a - 1. ) )**(-3.) if abs(a) > 1e-06 else 0.,
        a = 0., 
        b = a[mask], 
    )
    retval[mask] = integ * Ez[mask]

    if nu == 0:
        return retval

    # -- Linear growth rate calculations:     
    
    lnEfunc = lambda lna: -np.log( cosmo.efunc( np.exp(-lna) - 1 ) )
    retval  = (
        1. / ( retval * Ez**2 * a**2 ) 

        # 1-st log derivative of the hubble parameter function - this is calculated using a 4-point
        # finite difference method:   
        - reduce(
            lambda __sum, __w: ( 
                __sum + __w[0] * lnEfunc( np.log(a) + __w[1] * stepsize ) 
            ), 
            [(-1.0,  2.0), ( 8.0,  1.0), (-8.0, -1.0), ( 1.0, -2.0)], # finite difference weights and steps 
            np.zeros_like(a), 
        ) / ( 12*stepsize )
    )
    return retval

@dataclass
class _Settings:
    windowfunc: Literal["tophat", "gauss"] = "tophat"      
    lnka      : float                      = np.log(1e-06)
    lnkb      : float                      = np.log(1e+06)
    stepsize  : float                      = 0.01         

@dataclass(frozen = True)
class PowerSpectrum( ABC ):
    r"""
    Base class for matter power spectrum related calculations. This can be used for calculating matter 
    power spectrum and spectral moments.

    Note: Do not use this class directly, use the sub-classes. 

    Parameters
    ----------
    cosmo : astropy.cosmology.FLRW
        Cosmology model. 

    redshift : float
        Redshift at which the power spectrum is evaluated. This value must be greater than -1.

    ns : float, default=1
        Index of the initial power spectrum. This is assumed to be a power law function :math:`Ak^n_s`.

    sigma8: float, default=1
        Normalization of the power spectrum. Specify the value of the matter variance, smoothed by a 
        8 Mpc/h radius window.

    """
    
    cosmo   : FLRW
    redshift: float
    ns      : float     = 1.
    sigma8  : float     = 1.
    settings: _Settings = field( default_factory = _Settings, repr = False )
    
    # Other fields:
    norm    : float = field(default = 1.    , init = False, compare = False, repr = False) 
    dplus_z : float = field(default = np.nan, init = False, compare = False, repr = False) 
    dplus_0 : float = field(default = np.nan, init = False, compare = False, repr = False) 

    def __post_init__(self) -> None:
        object.__setattr__( self, 'dplus_0', linearGrowth(z = 0.           , cosmo = self.cosmo) )
        object.__setattr__( self, 'dplus_z', linearGrowth(z = self.redshift, cosmo = self.cosmo) )
        self.normalize()
        return
    
    @abstractmethod
    def transfer(
            self, 
            lnk: _T, 
        ) -> _T:
        r"""
        Calculate the value of the linear matter transfer function at current redshift.

        Parameters
        ----------
        lnk : array_like
            Scale at which the spectrum is evaluated. Natural logarithm of the wavenumber in 1/Mpc.

        Returns
        -------
        tf : array_like
            Transfer function values.  

        """
        raise NotImplementedError("transfer function model is not implemented")
    
    def power(
            self, 
            lnk: _T, 
            normalize: bool = True, 
        ) -> _T:
        r"""
        Calculate the value of the linear matter power spectrum at current redshift.

        Parameters
        ----------
        lnk : array_like
            Scale at which the spectrum is evaluated. Natural logarithm of the wavenumber in 1/Mpc.

        normalize : bool, default=True
            If true, normalize the value using sigma-8 parameter.

        Returns
        -------
        ps : array_like
            Power spectrum in Mpc^3 unit.  

        """

        # Initial power law spectrum:
        retval = np.exp( np.multiply(self.ns, lnk) ) * self.norm 
        
        # Power spectrum interpolated to current redshift:
        retval = retval * self.transfer(lnk)**2 

        if normalize:
            # Normalizing the power spectrum using the sigma-8 parameter value:
            retval = self.sigma8**2 * retval
        return retval
    
    def spectralMoment(
            self, 
            lnr       : _T, 
            nu        : Literal[0, 1] = 0, 
            j         : int           = 0,
            normalize : bool          = True,
        ) -> _T:
        r"""
        Calculate the j-th spectral moment or its log-derivative. 0-th spectral moment is the variance of
        matter density fluctuations (smoothed with the specified window function).

        Parameters
        ----------
        lnr : array_like
            Scale at which the value is evaluated. Natural logarithm of the radius in Mpc.

        nu : {0, 1}, optional
            Return the value of first log-derivative for non-zero value.
        
        j : int, default=0
            Order of the moment. Must be a non-negative integer.

        normalize : bool, default=True
            If true, normalize the value using sigma-8 parameter.

        Returns
        -------
        sj : array_like
            Spectral moment value or its first log derivative.  

        """

        def integrand(lnk: float, lnr: float, j: int, windowfunc: str) -> float:
            # Integrand for the spectral moment calculations: Delta(k) * Window(k*r)**2, where ``Delta(k)``
            # is the dimenssionless power spectrum and ``Window(x)`` is the smoothing function used to 
            # smooth out density fluctuations. 

            retval = self.power(lnk, normalize = False) * np.exp( np.multiply(3 + 2*j, lnk) )

            # Window functions: 
            kr = np.exp( np.add(lnk, lnr) )
            if windowfunc == "gauss":
                # Gaussian window function:
                retval *= np.exp( -kr**2 )
            else: # windowfunc == "tophat"
                # Spherical top-hat window function: ``3*j1(x) / x**2`` 
                retval *= ( 3*( np.sin(kr) - kr * np.cos(kr) ) / kr**3 )**2

            return retval
        
        def integrate(lnr: float, j: int, windowfunc: str) -> float:
            # Calculate the spectral moment integral for one radius value.
            retval, abserr = quad(
                integrand, 
                a    = self.settings.lnka, 
                b    = self.settings.lnkb, 
                args = ( lnr, j, windowfunc ),  
            )
            return retval if abs(retval) > 1e-08 else 0.
        
        integrateVectorized = np.vectorize(integrate, otypes = [np.float64])
        
        # -- Spectral moment calculation: 
        if nu == 0:
            retval = integrateVectorized(lnr, j, self.settings.windowfunc)
            if normalize: # using the sigma-8 parameter value...
                retval = self.sigma8**2 * retval
            return retval

        # -- Derivative calculation: 1-st log derivative of the spectral moment is calculated using a 
        # 4-point finite difference method:   
        lnsFunc = lambda lnr: np.log( integrateVectorized(lnr, j, self.settings.windowfunc) )
        retval  = (
            reduce(
                lambda __sum, __w: ( 
                    __sum + __w[0] * lnsFunc( lnr + __w[1] * self.settings.stepsize ) 
                ), 
                [(-1.0,  2.0), ( 8.0,  1.0), (-8.0, -1.0), ( 1.0, -2.0)], # finite difference weights and steps 
                np.zeros_like(lnr), 
            ) / ( 12*self.settings.stepsize )
        )
        return retval
    
    def normalize(self) -> None:
        r"""
        Normalise the power spectrum using the value of sigma-8 parameter.

        """
        # Save the current redshift and reset to present redshift (0), since the normalization
        # is calculated for present redshift... 
        redshift = self.redshift
        object.__setattr__(self, "redshift", 0.0)

        norm = 1. / self.spectralMoment(
            lnr        = np.log(8.) - np.log(self.cosmo.h), 
            nu         = 0, 
            j          = 0, 
            normalize  = False, 
        )
        object.__setattr__(self, "redshift", redshift) # reset redshift value
        return object.__setattr__( self, "norm", norm )
        
############################################################################################################
#                                   SOME MATTER POWER SPECTRUM MODELS!...                                  #
############################################################################################################

@dataclass(init = True, frozen = True, repr = False)
class EisensteinHu98_zb(PowerSpectrum):
    r"""
    A class implementing the power spectrum based on the linear matter transfer function by Eisenetein & Hu 
    (1998) with trace baryon content (no baryon oscillations).  

    Parameters
    ----------
    cosmo : astropy.cosmology.FLRW
        Cosmology model. 

    redshift : float
        Redshift at which the power spectrum is evaluated. This value must be greater than -1.

    ns : float, default=1
        Index of the initial power spectrum. This is assumed to be a power law function :math:`Ak^n_s`.

    sigma8: float, default=1
        Normalization of the power spectrum. Specify the value of the matter variance, smoothed by a 
        8 Mpc/h radius window.


    Usage
    -----
    >>> from astropy.cosmology import FlatLambdaCDM
    >>> cm = FlatLambdaCDM(H0 = 70., Om0 = 0.3, Ob0 = 0.05, Tcmb0 = 2.725)
    >>> ps = EisensteinHu98_zb(cm, redshift = 0., ns = 1., sigma8 = 1.)

    """

    # Internal paremeters:
    h      : float = field(init = False, default = np.nan, repr = False)
    Omh2   : float = field(init = False, default = np.nan, repr = False)
    theta  : float = field(init = False, default = np.nan, repr = False)
    z_eq   : float = field(init = False, default = np.nan, repr = False)
    z_d    : float = field(init = False, default = np.nan, repr = False)
    s      : float = field(init = False, default = np.nan, repr = False)
    alpha_g: float = field(init = False, default = np.nan, repr = False)

    def __post_init__(self) -> None:

        h    = self.cosmo.h          # Hubble parameter (in unit of 100 km/sec/Mpc)
        Omh2 = self.cosmo.Om0 * h**2 # Total matter density parameter 
        Obh2 = self.cosmo.Ob0 * h**2 # Density parameter for baryonic matter  
        object.__setattr__( self, 'h'   , h    )   
        object.__setattr__( self, 'Omh2', Omh2 )

        # CMB temperature in 2.7 K unit (if Tcmb0 is 0, default value of 2.725 is used). If using a astropy 
        # cosmology model or any model which stores this value as a Quantity object with a unit, np.array() 
        # will remove the the unit and get the value only: 
        theta = ( np.array(self.cosmo.Tcmb0) or 2.725 ) / 2.7
        object.__setattr__( self, 'theta', theta )

        # Redshift at matter-radiation equality (eqn. 1)
        z_eq = 2.5e+04 * Omh2 / theta**4
        object.__setattr__( self, 'z_eq', z_eq )

        # Redshift at drag epoch (eqn. 2)
        c1  = 0.313 * (1 + 0.607*Omh2**0.674) / Omh2**0.419
        c2  = 0.238*Omh2**0.223
        z_d = 1291.0*(Omh2**0.251)*(1 + c1*Obh2**c2) / (1 + 0.659*Omh2**0.828)
        object.__setattr__( self, 'z_d', z_d )

        # Sound horizon (eqn. 26)
        s = 44.5*np.log( 9.83/Omh2 ) / np.sqrt( 1 + 10.*Obh2**0.75 )
        object.__setattr__( self, 's', s )

        # Parameter alpha_Gamma, Eqn. 31
        fb = Obh2 / Omh2
        alpha_g = 1. - 0.328*np.log( 431*Omh2 ) * fb + 0.38*np.log( 22.3*Omh2 ) * fb**2
        object.__setattr__( self, 'alpha_g', alpha_g )
        
        return super().__post_init__()
    
    def transfer(
            self, 
            lnk: _T, /
        ) -> _T:
        Omh2, theta, alpha_g, s = self.Omh2, self.theta, self.alpha_g, self.s
        
        k     = np.exp(lnk) # wave-number in Mpc^-1
        gamma = Omh2 * ( alpha_g + ( 1 - alpha_g ) / ( 1 + ( 0.43*k*s )**4 ) ) # shape (Eqn 30)
        q     = k * ( theta**2 / gamma ) # dimension-less wavenumber 
        
        # Transfer function
        Lq = np.log( 2*np.exp(1) + 1.8*q )
        Cq = 14.2 + 731.0 / ( 1 + 62.5*q )
        retval = Lq / (Lq + Cq*q**2)
        
        # Interpolation to current redshift using linear growth factor 
        dplus_z = self.dplus_z / self.dplus_0 # linear growth factor, normalized to be 1 at redshift 0
        retval = dplus_z * retval 
        return retval


# Loook-up table for built-in models
availableModels: dict[str, type[PowerSpectrum]] = {
    "eisenstein98_zb": EisensteinHu98_zb, 
}


