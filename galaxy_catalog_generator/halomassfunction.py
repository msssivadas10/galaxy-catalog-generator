import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline
from typing import TypeVar, ClassVar, Literal, Any
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from .powerspectrum import PowerSpectrum

_T = TypeVar('_T')

@dataclass(frozen = True)
class HaloMassFunction( ABC ):
    r"""
    Base class for halo mass-function models.

    Parameters
    ----------
    psmodel : PowerSpectrum
        Matter power spectrum object. This is used for accessing the cosmology model as well as calculating
        the denisty variance values. It should be configured to the same redshift.
    
    redshift : float
        Redshift at which the mass-function is evaluated. This value must be greater than -1.

    Delta : float, default=200
        Halo over-density w.r.to the mean matter density. 

    """
    psmodel : PowerSpectrum
    redshift: float
    Delta   : int   = 200

    _sigtable: CubicSpline = field(default = None, init = False, repr = False)

    # Model id is used to used to identify the model and used as key for the model dict
    id: ClassVar[str] = ...  

    def __post_init__(self) -> None:
        assert isinstance(self.psmodel, PowerSpectrum)
        if not np.allclose(self.psmodel.redshift, self.redshift):
            raise ValueError( f"redshift mismatch: {self.psmodel.redshift} and {self.redshift}" )
        
        self.createInterpolationTable()
        return

    @property # Over-density threshold for spherical collapse
    def delta_sc(self) -> float: return 1.6864701998411453 

    @property
    def rho_m(self) -> float: 
        # Matter density at redshift 0
        
        import astropy.units as u

        rho_m =  u.Quantity(
            self.psmodel.cosmo.critical_density0 * self.psmodel.cosmo.Om0, 
            u.Msun / u.Mpc**3
        )
        return np.array(rho_m)

    def setRedshift(self, z: float) -> None:
        r"""
        Set the value of redshift.
        """
        self.psmodel.setRedshift(z) # set redshift power spectrum model
        object.__setattr__(self, "redshift", z)
        return self.createInterpolationTable() # recalculate tables

    @abstractmethod
    def fsigma(
            self, 
            s: _T, 
            **kwargs: Any, 
        ) -> _T:
        r"""
        Calcuate the halo mass-function value, given the variance.

        Parameters
        ----------
        s : array_like
            Mass variable - variance value corresponding the halo mass.

        **kwargs : Any
            Model specific keyword arguments.

        Returns
        -------
        fs : array_like
            Mass-function values.

        """
        pass

    def sigma(
            self, 
            lnm: _T, 
            return_derivative: bool = False, 
            interpolated: bool      = True,
        ) -> _T | tuple[_T, _T]:
        r"""
        Calculate the value of the matter variance or its first log derivative.

        Parameters
        ----------
        lnm : array_like
            Natural logarithm of the mass of the halo in Msun.

        return_derivative : bool, default=False
            If true, return the log derivative w.r.to mass also.
        
        interpolated : bool, default=True
            If true, calculate the values froma pre-calculated interpolation table. Otherwise, evaulate it 
            using the integral expression.

        Returns
        -------
        s : array_like
            Matter variance values.

        dlnsdlnm : array_like
            Log derivative - present only if ``return_derivative`` is true.

        """

        if interpolated and self._sigtable is not None:
            sigma = np.exp( self._sigtable(lnm) )
            if return_derivative:
                return sigma, self._sigtable(lnm, nu = 1)
            return sigma
        
        rho_m = self.rho_m # Matter density at redshift 0ated:
        rho_h = rho_m      # Halo density (TODO: chek if the halo density is rho_m * self.Delta)
        
        # Lagrangian radius corresponding to the mass
        lnr = ( lnm + np.log(3. / (4.*np.pi) / rho_h ) ) / 3.

        args = dict( j = 0, normalize = True, )

        sigma = np.sqrt( self.psmodel.spectralMoment( lnr, nu = 0, **args ) )
        if not return_derivative:
            return sigma
        
        dlnslnm = self.psmodel.spectralMoment( lnr, nu = 1, **args ) / 6.
        return sigma, dlnslnm
    
    def createInterpolationTable(
            self, 
            lnma   : float = np.log(1e+06), 
            lnmb   : float = np.log(1e+16), 
            tabsize: int   = 51, 
        ) -> None:
        r"""
        Calculate variance values and create an interpolation table.

        Parameters
        ----------
        lnma : float, default=log(1e+6)
            Minimum mass value for the interpolation table.

        lnmb : float, default=log(1e+16)
            Maximum mass value for the interpolation table.
        
        tabsize  : int, default=51
            Size of the initerpolation table.

        """

        lnm   = np.linspace( lnma, lnmb, tabsize )
        sigma = self.sigma( 
            lnm, 
            return_derivative = False, 
            interpolated      = False, 
        )
        sigtable = CubicSpline( lnm, np.log(sigma) )
        object.__setattr__( self, "_sigtable", sigtable )
        return

    def massFunction(
            self, 
            lnm: _T, 
            return_value: Literal["dndlogm", "dndlnm", "dndm", "fs"] = "dndlogm", 
        ) -> _T:
        r"""
        Calculate the halo mass-function in specified format.

        Parameters
        ----------
        lnm : array_like
            Natural logarithm of the mass of the halo in Msun.

        return_value : {'dndlogm', 'dndlnm', 'dndm', 'fs'}, optional
            Specify the format of the mass-function return value.

            'dndlogm'
                Number density of halos in unit base-10 logarithmic mass interval. Unit is ``1/Mpc^3``.
            'dndlnm'
                Number density of halos in unit logarithmic mass interval. Unit is ``1/Mpc^3``.
            'dndm'
                Number density of halos in unit mass interval. Unit is ``1/Mpc^3/Msun``.
            'fs'
                Dimenssionless mass-function. No units for this one.

        Returns
        -------
        hmf : array_like
            Halo mass-function.

        """

        s, dlnsdlnm = self.sigma(lnm, return_derivative = True)
        
        # Dimensionless mass-function, f(sigma): it does not have any unit
        retval = self.fsigma(s)
        if return_value == "fs":
            return retval
        
        m = np.exp(lnm) # mass in Msun

        # dn/dlnm: number density per unit log mass interval, unit:  1/Mpc3
        retval = retval * np.abs(dlnsdlnm) * self.rho_m / m
        if return_value == "dndlogm":
            return retval / 2.3025850929940455 
        elif return_value == "dndlnm":
            return retval 
        
        # dn/dm: number density per unit mass interval, unit: 1/Msun/Mpc3
        retval = retval / m
        return retval

############################################################################################################
#                                   LOAD FROM HALO MASS-FUNCTION TABLE                                     #
############################################################################################################
 
@dataclass(init = False, frozen = True, repr = False)
class MassFunctionData(HaloMassFunction):
    r"""
    Halo mass-function loaded from pre-calculated data. This can be given as an array or the path to the 
    text file containing the data (see ``data`` parameter). 

    Parameters
    ----------
    data : str, ndarray of shape (N, 2)
        Mass-function table as an array or path to the file containing the table. This should have two 
        columns - natural log of halo mass (Msun) and natural log of mass-function dn/dm (Mpc^-3).  
    
    psmodel : PowerSpectrum
        Matter power spectrum object. This is used for accessing the cosmology model as well as calculating
        the denisty variance values. It should be configured to the same redshift.
        
    redshift : float
        Redshift at which the mass-function is evaluated. This value must be greater than -1.

    Delta : float, default=200
        Halo over-density w.r.to the mean matter density. 

    **loaderargs : optional
        Other keyword arguments are passed to ``numpy.loadtxt`` for loading the data from file.

    """
    data: CubicSpline
    id: str 

    def __init__(
            self, 
            data    : NDArray[np.float64], 
            psmodel : PowerSpectrum,
            redshift: float,
            Delta   : int = 200,
        ) -> None:
        # ``data`` array should contain the natural log of halo mass (Msun) and natural log of 
        # mass-function dn/dm (Mpc^-3) as the first two columns.
        data = np.asarray(data)
        assert data.ndim == 2 and data.shape[1] >= 2
        x, y = data[:, 0], data[:, 1]    
        object.__setattr__(self, "data", CubicSpline(x, y)  )
        object.__setattr__(self, "id"  , f"data{ id(data) }")
        return super().__init__(psmodel = psmodel, redshift = redshift, Delta = Delta)
    
    def massFunction(
            self, 
            lnm: _T, 
            return_value: Literal["dndlogm", "dndlnm", "dndm", "fs"] = "dndlogm", 
        ) -> _T:
        # It is easier to directly calculate the mass-functions dn/dm or dn/dlnm than using f(s), which 
        # needed to be calculated from these values 

        lnm = np.asarray(lnm, dtype = np.float64)
        
        # Set the value outside the data range to 0;
        validrange         = np.logical_not( ( lnm < self.data.x[ 0] ) | ( lnm > self.data.x[-1] ) )
        retval             = np.zeros_like(lnm)
        retval[validrange] = np.exp( self.data(lnm[validrange]) )

        if return_value == "fs":
            # Calculate ``f(s)``, given mass and dn/dm.
            s, dlnsdlnm = self.sigma(lnm, return_derivative = True)
            retval      = retval * np.exp(2*lnm) / np.abs(dlnsdlnm) / self.rho_m
            return retval
        
        elif return_value != "dndm":
            retval = np.exp(lnm) * retval
            if return_value == "dndlnm" : return retval
            if return_value == "dndlogm": return retval * 2.3025850929940455 

        return retval
    
    def fsigma(
            self, 
            s: _T, 
        ) -> _T:
        # Calculating mass from given varience values, by finding roots of the spline:
        lnm  = self._sigtable.solve(s)
        return self.massFunction(lnm, return_value = "fs") 

############################################################################################################
#                                   SOME HALO MASS-FUNCTION MODELS!...                                     #
############################################################################################################

@dataclass(init = True, frozen = True, repr = False)
class MassFunctionTinker08(HaloMassFunction):
    r"""
    Halo mass-function by Tinker et al (2008). This model is cosmology independent, but has redshift and 
    over-density dependence.

    Parameters
    ----------
    psmodel : PowerSpectrum
        Matter power spectrum object. This is used for accessing the cosmology model as well as calculating
        the denisty variance values. It should be configured to the same redshift.

    redshift : float
        Redshift at which the mass-function is evaluated. This value must be greater than -1.

    Delta : float, default=200
        Halo over-density w.r.to the mean matter density. 

    """

    id: ClassVar[str] = "tinker08"

    # Internal parameters:
    A    : float = field(init = False, default = np.nan, repr = False)
    a    : float = field(init = False, default = np.nan, repr = False)
    b    : float = field(init = False, default = np.nan, repr = False)
    c    : float = field(init = False, default = np.nan, repr = False)
    alpha: float = field(init = False, default = np.nan, repr = False)

    # Table of parameters as function over-density, Delta in Tinker (2008)
    _params: ClassVar[CubicSpline] = CubicSpline(
        [   200,   300,   400,   600,   800,  1200,  1600,  2400,  3200 ] , # Delta   
        [[0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260 ] , # A
         [1.470, 1.520, 1.560, 1.610, 1.870, 2.130, 2.300, 2.530, 2.660 ] , # a
         [2.570, 2.250, 2.050, 1.870, 1.590, 1.510, 1.460, 1.440, 1.410 ] , # b
         [1.190, 1.270, 1.340, 1.450, 1.580, 1.800, 1.970, 2.240, 2.440 ]], # c
        axis = 1,
    )

    def __post_init__(self) -> None:
        A0, a0, b0, c = self._params(self.Delta)
        object.__setattr__( self, "c", c )

        alpha = pow(10.0, -( 0.75 / np.log10( self.Delta / 75.0 ) )**1.2) # Eqn. 8
        object.__setattr__( self, "alpha", alpha )
    
        zp1 = self.redshift + 1.
        A   = A0 * zp1**(-0.14 )     
        object.__setattr__( self, "A", A )

        a   = a0 * zp1**(-0.06 )     
        object.__setattr__( self, "a", a )

        b   = b0 * zp1**(-alpha)
        object.__setattr__( self, "b", b )

        return super().__post_init__()
    
    def fsigma(
            self, 
            s: _T, 
        ) -> _T:
        A, a, b, c = self.A, self.a, self.b, self.c        
        s = np.asarray(s, dtype = np.float64)
        f = A * (1 + (b / s)**a) * np.exp(-c / s**2)
        return f
    
availableModels: dict[str, type[HaloMassFunction]] = {
    "tinker08": MassFunctionTinker08,
}

