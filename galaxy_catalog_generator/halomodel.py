import numpy as np
from numpy.random import default_rng, Generator
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline
from scipy.special import erf
from typing import TypeVar, Literal
from dataclasses import dataclass, field
from astropy.cosmology import FLRW
from powerspectrum import PowerSpectrum, availableModels as powerspectrum_models
from halomassfunction import HaloMassFunction, availableModels as massfunction_models

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

    rng : Generator
        Random number generator to use.

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

    _cmtable: CubicSpline = field(default = None, init = False, repr = False)

    @classmethod
    def create(
            cls,
            Mmin     : float,
            sigmaM   : float,
            M0       : float,
            M1       : float,
            alpha    : float,
            scaleSHMF: float,
            slopeSHMF: float,
            redshift : float, 
            cosmo    : FLRW, 
            psmodel  : type[PowerSpectrum] | str    = "eisenstein98_zb", 
            mfmodel  : type[HaloMassFunction] | str = "tinker08",
            ns       : float = 1., 
            sigma8   : float = 1.,
            Delta    : int   = 200,
            seed     : int   = None, 
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
            Halo mass-function model to use.

        ns : float, default=1
            Index of the initial power spectrum. This is assumed to be a power law function :math:`Ak^n_s`.

        sigma8 : float, default=1
            Normalization of the power spectrum. Specify the value of the matter variance, smoothed by a 
            8 Mpc/h radius window.

        Delta : int, default=200
            Halo over-density w.r.to the mean matter density. 

        seed : int, optional 
            Seed value for the random number generator used.

        """
        
        assert psmodel in powerspectrum_models, f"unknown power spectrum model: {psmodel}"
        assert mfmodel in massfunction_models , f"unknown mass-function model: {mfmodel}"

        psmodel = powerspectrum_models.get(psmodel)(
            cosmo, 
            redshift = redshift, 
            ns       = ns, 
            sigma8   = sigma8, 
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
            psmodel   = psmodel, 
            mfmodel   = massfunction_models.get(mfmodel)(
                psmodel  = psmodel, 
                redshift = redshift, 
                Delta    = Delta, 
            ), 
            Delta    = Delta,
            rng      = default_rng(seed)
        )
        return self
    
    @property
    def delta_sc(self) -> float: return self.mfmodel.delta_sc

    @property 
    def growthFactor(self) -> float: 
        return self.psmodel.dplus_z / self.psmodel.dplus_0

    def createConcMassInterpolationTable(
            self, 
            crange : tuple[float, float] = (0., 2.), 
            tabsize: int = 64,
        ) -> None:
        r"""
        Calculate a table relating concentration `c` with the NFW mass function. This can be 
        used to find the inverse relation.

        Parameters
        ----------
        crange : (float, float), default=(0, 2)
            Specify the range of values for c to create the table.
        
        tabsize : int, default=64
            Size of the table.

        """
        c = np.linspace( crange[0], crange[1], tabsize )
        f = np.log(1 + c) - c / (1 + c)
        cmtable = CubicSpline(f, c)
        return object.__setattr__( self, "_cmtable", cmtable )
    
    def __post_init__(self) -> None:
        self.createConcMassInterpolationTable()
        return

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
        fsat[fsat < 0.] = 0.
        return ncen * fsat**self.alpha
    
    # TODO: subhaloMassFunction 
    
    def generateSatellitePositions(
            self, 
            lnm: float, 
        ) -> list[NDArray[np.float64]]:
        r"""
        Generate catalogs of satellite galaxies in a halo, given its mass. This catalog will have the X,Y,Z 
        position coordinates in Mpc (first 3 columns) and mass in Msun in the las column.

        Parameters
        ----------
        lnm : float
            Natural logarithm of the mass of the halo in Msun. 
            
            .. note::
                Sequence of values also accepted.

        Returns
        -------
        cat : ndarray
            Catalog of satellite galaxies in the halo. It can also be an empty array, if the halo do not 
            have any satellites. 

            .. note::
                If the input is a sequence of mass values, then return value will be a sequence of arrays, 
                corresponding the each halo.
        
        """

        lnm   = np.asarray(lnm, dtype = np.float64).flatten()
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

        haloMass   = np.exp(lnm)
        haloRadius = np.cbrt( 0.75 / np.pi * ( haloMass / rho_h ) ) # halo lagrangian radius in Mpc
        haloConc   = self.haloConcentration(lnm) # halo concentration parameter 

        positionList = []
        for Ns, radH, massH, concH in zip(satelliteN, haloRadius, haloMass, haloConc):
            if Ns <= 0: continue
            
            # Generating random values corresponding to the distance of the galaxy from the halo center.
            # This follows a distribution matching the NFW profile of the halo density. Sampling is done
            # using the inverse transformation method. 
            u     = self.rng.uniform(size = Ns)
            dist  = ( radH / concH ) * self._cmtable( u * ( np.log(concH + 1) - concH / (concH + 1) ) ) 
            theta = self.rng.uniform(size = Ns) *   np.pi
            phi   = self.rng.uniform(size = Ns) * 2*np.pi
            
            # Satellite galaxy coordinates x, y, and z in Mpc
            satellitePositions = np.array([
                dist * np.sin(theta) * np.cos(phi), 
                dist * np.sin(theta) * np.sin(phi),
                dist * np.cos(theta) 
            ]).T

            # Assigning random mass values to the satellite galaxies: These masses are drawn from a sub-
            # halo mass-function having the Schechter form (see <http://arxiv.org/abs/astro-ph/0402500v2>, 
            # Eqn 1). A rejection sampling is done to generate random values.
            massValues, remainingN = [], Ns
            while remainingN > 0:
                m = ( self.rng.pareto(self.slopeSHMF-1, size = remainingN) + 1. ) * self.Mmin / massH 
                # NOTE: acceptance probability is calculated for the above mentioned SHMF. If using a power
                # law with upper cut off, use the Pareto random number with proper scaling to match the 
                # wanted upper cutoff: that is, reject if mass is above the cutoff and no need to find
                # the acceptance probability...  
                massAccepted = m[ self.rng.uniform(size = remainingN) < np.exp( (1. - m) / self.scaleSHMF ) ]
                remainingN  -= len(massAccepted)
                massValues.append(massAccepted)
            massValues = np.concatenate(massValues) * massH

            positionList.append(
                np.concatenate(
                    [ satellitePositions, massValues[:,None] ], # mass in Msun along last column
                    axis = 1, 
                ) 
            )

        # positionList = np.concatenate(positionList, axis = 0)
        return positionList


