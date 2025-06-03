import logging
import numpy as np
from numpy.random import default_rng, Generator
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.special import erf
from typing import TypeVar, Literal
from dataclasses import dataclass, field
from astropy.cosmology import FLRW
from .powerspectrum import PowerSpectrum, availableModels as powerspectrum_models
from .halomassfunction import HaloMassFunction, MassFunctionData, availableModels as massfunction_models

_T = TypeVar('_T')

@dataclass(frozen = True)
class HaloModel:
    r"""
    A class for halo model calculations. This uses a 5-parameter halo occupation distribution model, with a 
    NFW density profile.

    Parameters
    ----------
    Mmin : float
        Specify the minimum mass (in Msun) for the halo to have at least one central galaxy. Same value is 
        used to check if the subhalo contain any (satellite) galaxy. 

    sigmaM : float
        Width of the central galaxy transition range. A value of 0 means the relation is a step function.

    M0 : float
        Minimum mass (in Msun) for halo to have satellite galaxies.
    
    M1 : float
        Scale factor for power law satelliete count relation (unit is Msun).

    alpha : float
        Index for the  power law satelliete count relation

    scaleSHMF : float
        Scale parameter for the subhalo mass-function. For a halo of mass M, ``m = beta*M`` gives the mass
        of the galaxy where the the function starts declining exponentially. A good range is [0.1, 0.5] and 
        ``beta < 1``.

    slopeSHMF : float
        Slope parameter for the subhalo mass-function. A value close to 2 can be used.

    redshift : float
        Redshift at which quantities are evaluated. Must be greater than -1.

    psmodel : PowerSpectrum
        Power spectrum model to use. Must be properly initialized.

    mfmodel : HaloMassFunction
        Halo mass-function model to use. Must be properly initialized.

    Delta : int, default=200
        Halo over-density w.r.to the mean matter density. 

    rng : Generator, optional
        Random number generator to use.

    meta : dict, optional
        Metadata related to the object. Can be used to store any values related to calculations.

    """
    Mmin     : float
    sigmaM   : float
    M0       : float
    M1       : float
    alpha    : float
    scaleSHMF: float
    slopeSHMF: float
    redshift : float 
    psmodel  : PowerSpectrum
    mfmodel  : HaloMassFunction
    Delta    : int = 200
    rng      : Generator = field(default_factory = default_rng, repr = False)
    meta     : dict      = field(default_factory = dict       , repr = False)

    _cmtable: CubicSpline = field(default = None, init = False, repr = False)

    @classmethod
    def create(
            cls,
            Mmin         : float,
            sigmaM       : float,
            M0           : float,
            M1           : float,
            alpha        : float,
            scaleSHMF    : float,
            slopeSHMF    : float,
            redshift     : float, 
            cosmo        : FLRW, 
            psmodel      : type[PowerSpectrum] | str    = "eisenstein98_zb", 
            mfmodel      : type[HaloMassFunction] | str = "tinker08",
            ns           : float = 1., 
            sigma8       : float = 1.,
            Delta        : int   = 200,
            seed         : int   = None,
            mf_loaderargs: dict  = None, 
            meta         : dict  = None,
        ):
        r"""
        Create a halo model object.

        Parameters
        ----------
        Mmin : float
            Specify the minimum mass (in Msun) for the halo to have at least one central galaxy.

        sigmaM : float
            Width of the central galaxy transition range. A value of 0 means the relation is a step function.

        M0 : float
            Minimum mass (in Msun) for halo to have satellite galaxies.
        
        M1 : float
            Scale factor for power law satelliete count relation (unit is Msun).

        alpha : float
            Index for the  power law satelliete count relation

        scaleSHMF : float
            Scale parameter for the subhalo mass-function. For a halo of mass M, ``m = beta*M`` gives the 
            mass of the galaxy where the the function starts declining exponentially. A good range is 
            [0.1, 0.5] and ``beta < 1``.

        slopeSHMF : float
            Slope parameter for the subhalo mass-function. A value close to 2 can be used.

        redshift : float
            Redshift.

        cosmo : FLRW
            Cosmology model to use.

        psmodel : subclass of ``PowerSpectrum``, str, default='eisenstein98_zb' 
            Power spectrum model to use.

        mfmodel : subclass of ``HaloMassFunction``, str, default='tinker08'
            Halo mass-function model to use. It can be the name of a pre-defined model or the path to the 
            file containing the values in plain text format 

        ns : float, default=1
            Index of the initial power spectrum. This is assumed to be a power law function :math:`Ak^n_s`.

        sigma8 : float, default=1
            Normalization of the power spectrum. Specify the value of the matter variance, smoothed by a 
            8 Mpc/h radius window.

        Delta : int, default=200
            Halo over-density w.r.to the mean matter density. 

        seed : int, optional 
            Seed value for the random number generator used.

        mf_loaderargs : dict
            Keyword arguments passed to ``numpy.loadtxt`` for loading the mass-function data from file.

        meta : dict, optional
            Metadata related to the object. Can be used to store any values related to calculations.

        """

        mf_loaderargs = mf_loaderargs or {}
        meta          = meta          or {}
        
        assert psmodel in powerspectrum_models, f"unknown power spectrum model: {psmodel}"
        powerspecObject = powerspectrum_models.get(psmodel)(
            cosmo, 
            redshift = redshift, 
            ns       = ns, 
            sigma8   = sigma8, 
        )
        
        if mfmodel in massfunction_models:
            # Create a mass-function object of this pre-defined model:
            massfuncObject = massfunction_models.get(mfmodel)(
                psmodel  = powerspecObject, 
                redshift = redshift, 
                Delta    = Delta, 
            )
        else:
            # This may be the filename containg the pre-calculated data or data...
            massfuncObject = MassFunctionData(
                mfmodel, 
                psmodel  = powerspecObject, 
                redshift = redshift, 
                Delta    = Delta, 
                **mf_loaderargs
            )

        self = HaloModel(
            Mmin      = Mmin,
            sigmaM    = sigmaM,
            M0        = M0,
            M1        = M1,
            alpha     = alpha,
            scaleSHMF = scaleSHMF, 
            slopeSHMF = slopeSHMF, 
            redshift  = redshift, 
            psmodel   = powerspecObject, 
            mfmodel   = massfuncObject, 
            Delta     = Delta,
            rng       = default_rng(seed),
            meta      = meta,
        )
        return self
    
    @property
    def delta_sc(self) -> float: return self.mfmodel.delta_sc

    @property 
    def growthFactor(self) -> float: 
        return self.psmodel.dplus_z / self.psmodel.dplus_0
    
    @property
    def cosmo(self) -> FLRW: return self.psmodel.cosmo
    
    def setRedshift(
            self, 
            z: float, 
        ) -> None:
        r"""
        Set the value of redshift. Value should be greater than -1.

        """
        # Set the redshift value to mass-function model: this will also set the redshift
        # of the power spectrum model. 
        self.mfmodel.setRedshift(z)
        return object.__setattr__(self, "redshift", z)
    
    def _updateHaloParameters(
            self, 
            **kwargs: float, 
        ) -> None:
        r"""
        Set new values for halo model parameters. Only values for the parameters ``Mmin``, ``sigmaM``, 
        ``M0``, ``M1``, ``alpha``, ``scaleSHMF`` and ``slopeSHMF`` can be updated.

        """
        updatableFields = ( "Mmin", "sigmaM", "M0", "M1", "alpha", "scaleSHMF", "slopeSHMF" )
        for __field, __value in kwargs.items():
            if __field not in updatableFields:
                raise AttributeError(f"cannot update value for attribute {__field!r}")
            object.__setattr__(self, __field, __value)
        return

    def createConcMassInterpolationTable(
            self, 
            crange : tuple[float, float] = (-3., 2), 
            tabsize: int = 64,
        ) -> None:
        r"""
        Calculate a table relating concentration `c` with the NFW mass function. This can be 
        used to find the inverse relation.

        Parameters
        ----------
        crange : (float, float), default=(0, 2)
            Specify the range of values for log10(c) to create the table.
        
        tabsize : int, default=64
            Size of the table.

        """
        c = np.logspace( crange[0], crange[1], tabsize )
        f = np.log(1 + c) - c / (1 + c)
        cmtable = CubicSpline(f, c)
        return object.__setattr__( self, "_cmtable", cmtable )
    
    def __post_init__(self) -> None:
        self.setRedshift( self.redshift )
        self.createConcMassInterpolationTable()
        return
    
    ########################################################################################################
    #                                               MODELS                                                 #
    ########################################################################################################

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
        return self.mfmodel.sigma(lnm, return_derivative = return_derivative, interpolated = interpolated)
    
    def haloMassFunction(
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
        return self.mfmodel.massFunction(lnm, return_value = return_value)
    
    def haloConcentration(
            self, 
            lnm: _T,     
        ) -> _T:
        r"""
        Calculate the concentration parameter for a halo of given mass. This value if calculated for the 
        currently set redshift value.

        Parameters
        ----------
        lnm : array_like
            Natural logarithm of the mass of the halo in Msun.

        Returns
        -------
        c : array_like
            Halo concentraion parameter.

        """
        # Redshift dependent parameters:
        a  = 1. / (1. + self.redshift)
        c0 = 3.395 * a**( 0.215 )
        b  = 0.307 * a**(-0.540 )
        g1 = 0.628 * a**( 0.047 )
        g2 = 0.317 * a**( 0.893 )
        v0 = ( 
                4.135 - 0.564 / a - 0.210 / a**2 + 0.0557 / a**3 - 0.00348 / a**4 
             ) / self.growthFactor
        v  = self.delta_sc / self.sigma(lnm, return_derivative = False)
        t  = v / v0

        # Concentration-mass relation:
        c  = c0 * t**(-g1) * (1. + t**(1./b))**(-b*(g2 - g1))
        return c
    
    def centralCount(
            self, 
            lnm: _T, 
        ) -> _T:
        r"""
        Return the average count of central galaxies in a halo, given its mass. This will be a sigmoid  
        function with smoothness controlled by the ``sigmaM`` parameter. If it is 0, then it will be a 
        step function.

        Parameters
        ----------
        lnm : array_like
            Natural logarithm of the mass of the halo in Msun.

        Returns
        -------
        Ncen : array_like
            Average number of central galaxies.
            
        """
        
        lnm = np.asarray(lnm, dtype = np.float64)
        if abs(self.sigmaM) < 1e-06:
            return np.heaviside( lnm - np.log(self.Mmin), 1. )
         
        x = ( lnm - np.log(self.Mmin) ) / self.sigmaM
        return 0.5 * ( 1. + erf(x) )

    def satelliteCount(
            self, 
            lnm: _T, 
        ) -> _T:
        r"""
        Return the average count of satellite galaxies in a halo, given its mass. This will be a power law
        specified by the parameters ``(M0, M1, alpha)``. 

        Parameters
        ----------
        lnm : array_like
            Natural logarithm of the mass of the halo in Msun.

        Returns
        -------
        Nsat : array_like
            Average number of satellite galaxies.
            
        """
        ncen = self.centralCount(lnm)
        fsat = (np.exp(lnm) - self.M0) / self.M1
        fsat = np.where(fsat < 0., 0., fsat)
        return ncen * fsat**self.alpha
    
    def subhaloMassFunction(
            self, 
            x: _T,
            lnm: float,
        ) -> _T:
        r"""
        Calculate the subhalo mass-function for given halo mass.

        Parameters
        ----------
        x : array_like
            Mass of the subhalo as fraction of the parent halo.

        lnm : float
            Natural logarithm of the mass of the halo in Msun. 

        Returns
        -------
        shmf : array_like
            Value of the subhalo mass-function, without normalization. This correspond to the probability 
            distribution of the sub-halo mass values. Multiply this with the correct normalization factor 
            to get the actual mass-function.

        """
        xMin, xMax = (self.Mmin / np.exp(lnm)), self.scaleSHMF
        slope      = self.slopeSHMF
        factor     = slope*xMin**slope / ( 1 - (xMin / xMax)**slope )

        x = np.asarray(x, dtype = np.float64)
        validRange = ( x >= xMin ) & ( x <= xMax )

        y = np.zeros_like(x)
        y[ validRange ] = factor * x[ validRange ]**(-slope - 1) 
        return y
    
    ########################################################################################################
    #                                       GALAXY CATALOG GENERATION                                      #
    ########################################################################################################
    
    def generateSatellitePositions(
            self, 
            lnm: float, 
            pos: tuple[float, float, float] = None,
        ) -> NDArray[np.float64]:
        r"""
        Generate catalogs of satellite galaxies in a halo, given its mass. This catalog will have the X,Y,Z 
        position coordinates in Mpc (first 3 columns) and mass in Msun in the las column.

        Parameters
        ----------
        lnm : float
            Natural logarithm of the mass of the halo in Msun. 

        pos : tuple[float, float, float], optional
            Position of the halo in Mpc. If given, the catalog also contains central galaxies, with position  
            and mass same as the halo.

        Returns
        -------
        cat : ndarray or None
            Catalog of satellite galaxies in the halo. It can also be ``None`` value, if the halo do not 
            have any satellites. If optional halo position is given, first item is the central galaxy and 
            the rest (if any) are satellites. If not, only satellites are returned and the position will 
            be relative to that of halo.  
        
        """

        logger = logging.getLogger()

        lnm   = np.asarray(lnm, dtype = np.float64)#.flatten()
        z     = self.redshift
        rho_m = self.mfmodel.rho_m * (1 + z)**3 # Matter density at redshift z
        rho_h = rho_m # Halo density (TODO: chek if the halo density is rho_m * self.Delta)

        # Central galaxy count: this is drawn from a binomial distribution of n = 1, so that the value is 
        # either 1 or 0. If the model uses a step function use the average count as the actual count. 
        centralN = self.centralCount(lnm)
        if abs( self.sigmaM ) > 1e-06:
            centralN = self.rng.binomial(
                1, 
                p    = centralN, 
                size = np.shape(lnm), 
            )

        # Satellite galaxy count: this drawn from a poisson distribution with the calculated average
        satelliteN = self.rng.poisson(
            lam  = self.satelliteCount(lnm), 
            size = np.shape(lnm),
        )
        logger.debug(f"using central galaxy count={centralN} and satellite count={satelliteN}")

        haloMass   = np.exp(lnm)
        haloRadius = np.cbrt( 0.75 / np.pi * ( haloMass / rho_h ) ) # halo lagrangian radius in Mpc
        haloConc   = self.haloConcentration(lnm) # halo concentration parameter 

        pos        = None if pos is None else np.asarray(pos)
        Nc, Ns, radH, massH, concH, posH = centralN, satelliteN, haloRadius, haloMass, haloConc, pos
        if Ns <= 0: 
            # If the halo has a central galaxy, add that:
            if Nc > 0. and posH is not None:
                return np.concatenate(
                        [ [posH], [[massH]] ], # mass in Msun along last column
                        axis = 1, 
                    )
            logger.debug(f"no satellites for halo mass={haloMass:.3e} Msun") 
            return np.empty(shape = (0, 4), dtype = np.float64)
        
        # Generating random values corresponding to the distance of the galaxy from the halo center.
        # This follows a distribution matching the NFW profile of the halo density. Sampling is done
        # using the inverse transformation method. 
        u     = self.rng.uniform(size = Ns)
        dist  = ( radH / concH ) * self._cmtable( u * ( np.log(concH + 1) - concH / (concH + 1) ) ) 
        theta = np.arccos( self.rng.uniform(-1., 1., size = Ns) )
        phi   = 2*np.pi * self.rng.uniform(size = Ns)
        
        # Satellite galaxy coordinates x, y, and z in Mpc
        satellitePositions = np.array([
            dist * np.sin(theta) * np.cos(phi), 
            dist * np.sin(theta) * np.sin(phi),
            dist * np.cos(theta) 
        ]).T
        if posH is not None:
            satellitePositions = posH + satellitePositions

        # Assigning random mass values to the satellite galaxies: These masses are drown from a bounded 
        # pareto distribution, with bounds [Mmin, scaleSHMF*massH] and slope given by slopeSHMF. 
        # NOTE: using inverse transform sampling <https://en.wikipedia.org/wiki/Pareto_distribution>
        if self.scaleSHMF <= self.Mmin / massH:
            # In this case, the mass range is non existent, meaning that there are no satellites.
            logger.debug(f"ignoring satellites - minimum mass fraction {self.Mmin / massH:.3g} >= {self.scaleSHMF:.3g}")
            satellitePositions = np.empty(shape = (0, 3), dtype = satellitePositions.dtype)
            massValues         = np.array([])
        else:
            Ha, La     = self.scaleSHMF**self.slopeSHMF, (self.Mmin / massH)**self.slopeSHMF
            massValues = massH * ( 
                -( self.rng.uniform(size = Ns) * ( Ha - La ) - Ha ) / ( Ha * La )
            )**( -1. / self.slopeSHMF )

        # Another choise for the mass distribution will be a subhalo mass-function having the Schechter  
        # form (see <http://arxiv.org/abs/astro-ph/0402500v2>, Eqn 1). This requires rejection sampling 
        # with a pareto distribution as the initial distribution.
        # massValues, remainingN = [], Ns
        # while remainingN > 0:
        #     m = ( self.rng.pareto(self.slopeSHMF-1, size = remainingN) + 1. ) * self.Mmin / massH 
        #     massAccepted = m[ self.rng.uniform(size = remainingN) < np.exp( (1. - m) / self.scaleSHMF ) ]
        #     remainingN  -= len(massAccepted)
        #     massValues.append(massAccepted)
        # massValues = np.concatenate(massValues) * massH

        # If halo position is also given, then add a central galaxy with same mass as the halo: this 
        # will be the first item in the list...
        if posH is not None:
            satellitePositions = np.vstack([ posH , satellitePositions ])
            massValues         = np.hstack([ massH, massValues ])

        galaxyData = np.concatenate(
            [ satellitePositions, massValues[:,None] ], # mass in Msun along last column
            axis = 1, 
        )
        logger.debug(f"galaxy generation complete: {galaxyData.shape[0]} galaxies generated!") 
        return galaxyData 
    
    ########################################################################################################
    #                                   AVERAGE DENSITY CALCULATION                                        #
    ########################################################################################################
    
    def averageHaloDensity(
            self, 
            lnma: float = None,
            lnmb: float = np.log(1e+18),
        ) -> float:
        r"""
        Return the average halo number density at current redshift.

        Parameters
        ----------
        lnma : float, optional
            Lower limit for halo mass (natural log of value in Msun). If not specified, a value 4 sigma 
            units lower than the minimum mass for galaxy formation is used.

        lnmb : float, default=log(1e+18)
            Upper limit for halo mass (natural log of value in Msun).

        Returns
        -------
        retval : float
            Halo number density in Mpc^-3.

        """

        # Integration limits: lower limit is calculated 4 sigma units below the minimum mass.
        a = lnma or np.log( self.Mmin ) - 4*self.sigmaM
        b = lnmb
        assert a < b, "lower limit must be less than upper limit"

        retval, abserr = quad(
            self.haloMassFunction, 
            a    = a, 
            b    = b,
            args = ( "dndlnm", )
        )
        return retval

    def averageGalaxyDensity(
            self, 
            lnma: float = None,
            lnmb: float = np.log(1e+18),
        ) -> float:
        r"""
        Return the average galaxy number density at current redshift.

        Parameters
        ----------
        lnma : float, optional
            Lower limit for halo mass (natural log of value in Msun). If not specified, a value 4 sigma 
            units lower than the minimum mass for galaxy formation is used.

        lnmb : float, default=log(1e+18)
            Upper limit for halo mass (natural log of value in Msun).

        Returns
        -------
        retval : float
            Galaxy number density in Mpc^-3.
             
        """

        # Integration limits: lower limit is calculated 4 sigma units below the minimum mass.
        a = lnma or np.log( self.Mmin ) - 4*self.sigmaM
        b = lnmb
        assert a < b, "lower limit must be less than upper limit"

        def integrand(lnm: float) -> float:
            galaxyCount = self.centralCount(lnm) + self.satelliteCount(lnm)
            return galaxyCount * self.haloMassFunction(lnm, return_value = "dndlnm")

        retval, abserr = quad( integrand, a = a, b = b )
        # retval /= self.averageHaloDensity(lnma = lnma, lnmb = lnmb) # normalizing with halo density
        return retval

    def averageSatelliteFraction(
            self, 
            lnma: float = None,
            lnmb: float = np.log(1e+18),
        ) -> float:
        r"""
        Return the average satellite galaxy number density at current redshift, as a fraction of the halo 
        number density.

        Parameters
        ----------
        lnma : float, optional
            Lower limit for halo mass (natural log of value in Msun). If not specified, a value 4 sigma 
            units lower than the minimum mass for galaxy formation is used.

        lnmb : float, default=log(1e+18)
            Upper limit for halo mass (natural log of value in Msun).

        Returns
        -------
        retval : float
            Satellite galaxy number fraction.
             
        """

        # Integration limits: lower limit is calculated 4 sigma units below the minimum mass.
        a = lnma or np.log( self.Mmin ) - 4*self.sigmaM
        b = lnmb
        assert a < b, "lower limit must be less than upper limit"

        def integrand(lnm: float) -> float:
            galaxyCount = self.satelliteCount(lnm)
            return galaxyCount * self.haloMassFunction(lnm, return_value = "dndlnm")

        retval, abserr = quad( integrand, a = a, b = b )
        # retval /= self.averageHaloDensity(lnma = lnma, lnmb = lnmb) # normalizing with halo density
        return retval / self.averageGalaxyDensity(lnma, lnmb)
    


