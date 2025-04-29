
__all__ = ["HODOptimizer"]

import logging
import numpy as np, numpy.typing as nt
from dataclasses import dataclass, field, asdict
from scipy.optimize import minimize
from ..halomodel import HaloModel

@dataclass(frozen = True)
class HODOptimizerResult:
    Mmin   : float
    M1     : float 
    galaxy_density    : float 
    satellite_fraction: float 
    score  : float 
    nit    : int 
    status : int 
    message: str 
    success: bool
    def asdict(self) -> dict: return asdict(self)     

@dataclass(frozen = True)
class HODOptimizer:
    r"""
    An object for optimizing the halo models. This will find the optimum value for the parameters ``Mmin`` 
    and ``M1`` for a given galaxy density and satellite fraction values. Parameter ``M0`` is set to the 
    same as ``Mmin``. 

    Parameters
    ----------
    halomodel : HaloModel
        Partially initialised halo model object. 

    galaxy_density : float 
        Required value for galaxy density at current redshift, in Mpc^-3 units.

    satellite_fraction : float
        Required value for satellite fraction at current redshift.

    Mmin_range : tuple[float, float], default=( 1e+12, 1e+14 )
        Search range for the parameter ``Mmin``. Value is in Msun units.
    
    M1_range : tuple[float, float], default=( 1e+13, 1e+15 )
        Search range for the parameter ``M1``. Value is in Msun units.

    gridsize : int, default=12
        Size of the grid used to get initial guess.

    Usage
    -----
    >>> import astropy.cosmology as cm
    >>> hm = HaloModel.create(Mmin = -1, sigmaM = 0., M0 = -1, M1 = -1, alpha = 1., scaleSHMF = 0.5,
    ...         slopeSHMF = 2.0, cosmo = cm.Planck18, redshift  = 3., psmodel = "eisenstein98_zb", 
    ...         mfmodel = "tinker08", ns = 0.9649, sigma8 = 0.807952, Delta = 200)
    >>> optimizer = HODOptimizer(hm, galaxy_density = 2.5e-4, satellite_fraction = 0.05, 
    ...                 Mmin_range = [ 1e+11, 1e+14 ], M1_range = [ 1e+13, 1e+15 ])
    >>> res = optimizer.optimize()
    >>> res.Mmin
    1151736838496.8633
    >>> res.M1
    22431431844829.793
    >>> res.galaxy_density
    0.00025000000260087286
    >>> res.satellite_fraction
    0.050000000096677676
    
    """
    halomodel         : HaloModel
    galaxy_density    : float 
    satellite_fraction: float
    Mmin_range        : tuple[float, float] = ( 1e+11, 1e+14 )
    M1_range          : tuple[float, float] = ( 1e+13, 1e+15 )
    gridsize          : int  = field( default = 12, repr = False )    
    _cache            : dict = field( init = False, repr = False, default_factory = dict )   

    def score(
            self,
            p: tuple[float, float], 
        ) -> float:
        r"""
        A estimate for the weighted difference of the calculated and required galaxy density and satellite 
        fraction values.

        Parameters
        ----------  
        p : tuple[float, float]
            Point on the ``ln(Mmin)`` - ``ln(M1)`` plane where the function is evaluated.

        Returns
        -------
        score : float
            Log of the weighted difference. 

        """

        # Set the HOD paramters to the model
        Mmin, M1  = np.exp(p)
        halomodel = self.halomodel
        halomodel._updateHaloParameters(Mmin = Mmin, M1 = M1, M0 = Mmin)

        # Relative diffrence in galaxy density
        galaxy_density = halomodel.averageGalaxyDensity( lnmb = np.log(1e+18) )
        deltaG         = ( galaxy_density / self.galaxy_density - 1. )**2
        
        # Relative difference in satellite fraction
        satellite_fraction = halomodel.averageSatelliteFraction( lnmb = np.log(1e+18) )
        deltaS             = ( satellite_fraction / self.satellite_fraction - 1. )**2
        
        # Total score: weighted distance from observed and calculated values
        score = np.log( deltaG + deltaS + 1e-16 ) 

        self._cache['galaxy_density']     = galaxy_density
        self._cache['satellite_fraction'] = satellite_fraction
        self._cache['score']              = score
        return score
    
    def optimize(
            self, 
            p0: tuple[float, float] = None,
        ) -> HODOptimizerResult:
        r"""
        Optimize the HOD parameters in the halo model.

        Parameters
        ----------  
        p0 : tuple[float, float], optional
            Initial guess for the minimum point.

        Returns
        -------
        res : HODOptimizerResult
            An object storing details of the minimization process.

        """
        logger = logging.getLogger()
        
        x0, y0 = p0 if p0 is not None else self.gridOptimize()
        xa, xb = self.Mmin_range
        ya, yb = self.M1_range
        
        # Minimizing the score function to get the optimum values: 
        logger.debug(f"optimizing for Mmin and M1 with guess ({x0=:.3f}, {y0=:.3f})...")
        logger.debug(f"using Mmin range [{xa:.3e}, {xb:.3e}] and M1 range [{ya:.3e}, {yb:.3e}]...")
        res = minimize(
            self.score, 
            x0     = [ x0, y0 ], 
            bounds = [
                ( np.log(xa), np.log(xb) ), # Mmin range
                ( np.log(ya), np.log(yb) ), # M1 range
            ], 
            args   = (),  
        )
        if not res.success:
            logger.warning( f"optimization failed with messagee {res.message!r} after {res.nit} iterations" )

        # Set the optimum values to the halo model
        score = self.score(res.x)
        logger.debug(f"optimum values: Mmin={self.halomodel.Mmin:.3e} Msun, M1={self.halomodel.M1:.3e} Msun")
        logger.debug(f"Final score: {score:.4g}")
        logger.debug(f"Galaxy density: {self._cache['galaxy_density']:.3e} Mpc^-3")
        logger.debug(f"Satellite fraction: {self._cache['satellite_fraction']:.3g}") 
        
        res = HODOptimizerResult(
            Mmin               = self.halomodel.Mmin, 
            M1                 = self.halomodel.M1, 
            galaxy_density     = self._cache["galaxy_density"], 
            satellite_fraction = self._cache["satellite_fraction"],
            score              = self._cache["score"], 
            nit                = res.nit, 
            status             = res.status, 
            message            = res.message,
            success            = res.success,
        )
        self._cache.clear()
        return res
    
    @property
    def xvalues(self) -> nt.NDArray[np.float64]: # log(Mmin) values
        xa, xb = self.Mmin_range
        return np.linspace( np.log( xa ), np.log( xb ), self.gridsize )
    
    @property
    def yvalues(self) -> nt.NDArray[np.float64]: # log(M1) values
        ya, yb = self.M1_range
        return np.linspace( np.log( ya ), np.log( yb ), self.gridsize )
    
    def gridOptimize(self) -> tuple[float, float]:
        r"""
        Minimize the score function using a grid. 
        
        This will calculate the score function on a grid in the Mmin - M1 plane, using the bounds and size 
        provided and get the point where the value is minimum. This can be used to get approximate values 
        of the parameters ``ln(Mmin)`` and ``ln(M1)``, to use as initial guess for a more accurate value 
        using minimization.

        Returns
        -------
        x : float
            Natural log of the parameter ``Mmin`` in Msun.

        y : float 
            Natural log of the parameter ``M1`` in Msun.

        """
        logger = logging.getLogger()

        # Creating the grid: this is used to locate the initial guess for the minimum
        logger.debug("generating grid...")
        x, y  = np.meshgrid(self.xvalues, self.yvalues)
        score = np.zeros_like(x)
        for i in range( score.shape[0] ):
            for j in range( score.shape[1] ):
                score[i, j] = self.score([ x[i, j], y[i, j] ])
        
        # Initial guess is the minimum point in the grid:
        logger.debug("finding minimum from grid...")
        i, j   = np.unravel_index( np.argmin(score), score.shape )
        x0, y0 = x[i, j], y[i, j]

        logger.debug(f"minimum from grid: ln(Mmin)={x0} and ln(M1)={y0}")
        return x0, y0
    