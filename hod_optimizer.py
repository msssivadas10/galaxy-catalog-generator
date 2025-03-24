__version__ = "0.1a"

import os, os.path
import numpy as np
from scipy.optimize import minimize
from astropy.cosmology import w0waCDM
from galaxy_catalog_generator.halomodel import HaloModel

# Look-up table to get `sigma_8`values used in the abacus summit cosmologies. This data is not in the
# catalog headers, so it is taken from <https://abacussummit.readthedocs.io/en/latest/cosmologies.html>.
_AbacusCosmologySigma8Table = {
    "000": 0.807952, "001": 0.776779, "002": 0.808189, "003": 0.855190, "004": 0.749999, "009": 0.811362, 
    "010": 0.823630, "011": 0.905993, "012": 0.827899, "013": 0.813715, "014": 0.800000, "015": 0.835005, 
    "016": 0.793693, "017": 0.815903, "018": 0.819708, "019": 0.805050, "020": 0.811350, "021": 0.849842,   
    "022": 0.824961, "100": 0.808181, "101": 0.808156, "102": 0.808270, "103": 0.808075, "104": 0.808166, 
    "105": 0.808177, "106": 0.808181, "107": 0.808168, "108": 0.808177, "109": 0.808161, "110": 0.808170, 
    "111": 0.808179, "112": 0.824120, "113": 0.792107, "114": 0.808245, "115": 0.808089, "116": 0.866133, 
    "117": 0.808205, "118": 0.808152, "119": 0.808168, "120": 0.808181, "121": 0.808169, "122": 0.808170,   
    "123": 0.808176, "124": 0.808175, "125": 0.811995, "126": 0.803908, "130": 0.711201, "131": 0.807866, 
    "132": 0.808189, "133": 0.905163, "134": 0.808458, "135": 0.808160, "136": 0.882475, "137": 0.729693, 
    "138": 0.821432, "139": 0.793003, "140": 0.772033, "141": 0.799575, "142": 0.707082, "143": 0.847522, 
    "144": 0.891360, "145": 0.801404, "146": 0.762696, "147": 0.777157, "148": 0.714913, "149": 0.824854, 
    "150": 0.937655, "151": 0.860820, "152": 0.677885, "153": 0.794389, "154": 0.838698, "155": 0.735302, 
    "156": 0.801974, "157": 0.872315, "158": 0.829816, "159": 0.718521, "160": 0.876756, "161": 0.793066, 
    "162": 0.779589, "163": 0.838824, "164": 0.774159, "165": 0.835954, "166": 0.837463, "167": 0.768419, 
    "168": 0.871407, "169": 0.777925, "170": 0.716059, "171": 0.852878, "172": 0.765650, "173": 0.763962, 
    "174": 0.840113, "175": 0.708760, "176": 0.892483, "177": 0.806026, "178": 0.791239, "179": 0.775969, 
    "180": 0.894071, "181": 0.730036, 
}

def _buildFromAbacusSimulation(
        simname: str, 
        redshift: float,
        M0           : float = None, 
        sigmaM       : float =  0.0, 
        alpha        : float =  1.0, 
        scaleSHMF    : float =  0.5, 
        slopeSHMF    : float = -2.0, 
        powerspectrum: str  = "eisenstein98_zb", 
        massfunction : str  = "tinker08",  
    ) -> HaloModel:

    import re
    from abacusnbody.metadata import get_meta

    # Get the simulation parameters using simulation name:
    header = get_meta(simname, redshift)
    header["CMBTemperature"] = 2.728

    # Sigma-8 value is taken from the look-up table, using the cosmology ID
    header["sigma8_m"] = _AbacusCosmologySigma8Table.get(
            re.search(
                r"AbacusSummit_\w+_c(\d\d\d)_ph\d\d\d",  
                header["SimName"], 
            ).group(1), 
            None, 
        )
    
    # Building the initial model, without HOD parameters Mmin and M1...
    halomodel = HaloModel.create(
        Mmin      = -1,       # Mmin is calculated later...
        sigmaM    = sigmaM, 
        M0        = M0 or -1, # Use M0 = Mmin, if not specified
        M1        = -1,       # M1 is calculated later...
        alpha     = alpha, 
        scaleSHMF = scaleSHMF, 
        slopeSHMF = slopeSHMF, 
        redshift  = header["Redshift"], 
        cosmo     = w0waCDM(
            H0    = header["H0"], 
            Om0   = header["Omega_M"], 
            Ode0  = header["Omega_DE"],
            Ob0   = header["omega_b"] / ( 0.01*header["H0"] )**2,
            w0    = header["w0"], 
            wa    = header["wa"], 
            Tcmb0 = header["CMBTemperature"],
        ), 
        psmodel = powerspectrum, 
        mfmodel = massfunction, 
        ns      = header["n_s"], 
        sigma8  = header["sigma8_m"], 
        Delta   = header["SODensity"][0],
    )
    return halomodel

def _scoreFunction(
        x: tuple[float, float], 
        halomodel         : HaloModel,
        galaxy_density    : float, 
        satellite_fraction: float,
        return_model      : bool = False, 
    ) -> float:

    # Set the HOD paramters t the model
    Mmin, M1 = np.exp(x)
    object.__setattr__(halomodel, "Mmin", Mmin)
    object.__setattr__(halomodel, "M1"  , M1  )
    if halomodel.M0 < 0: 
        object.__setattr__(halomodel, "M0"  , Mmin)

    # Relative diffrence in galaxy density
    deltaG = ( halomodel.averageGalaxyDensity( lnmb = np.log(1e+18) ) / galaxy_density - 1. )**2
    
    # Relative difference in satellite fraction
    deltaS = ( halomodel.averageSatelliteFraction( lnmb = np.log(1e+18) ) / satellite_fraction - 1. )**2
    
    # Total score: weighted distance from observed and calculated values
    score = np.log( deltaG + deltaS + 1e-16 ) 
    if return_model:
        return halomodel
    return score

def optimizeHaloModel(
        halomodel         : HaloModel,
        galaxy_density    : float, 
        satellite_fraction: float,
        Mmin_range        : tuple[float, float] = ( 1e+12, 1e+14 ),
        M1_range          : tuple[float, float] = ( 1e+13, 1e+15 ),
        gridsize          : int = 12, 
    ) -> HaloModel:

    # Creating the grid: this is used to locate the initial guess for the minimum
    ( xa, xb ), ( ya, yb ) = Mmin_range, M1_range
    x, y = np.meshgrid(
        np.linspace( np.log( xa ), np.log( xb ), gridsize ), # Mmin values
        np.linspace( np.log( ya ), np.log( yb ), gridsize ), # M1 values
    )
    scoreGrid = np.zeros_like(x)
    for i in range( scoreGrid.shape[0] ):
        for j in range( scoreGrid.shape[1] ):
            scoreGrid[i, j] = _scoreFunction(
                ( x[i, j], y[i, j] ), 
                halomodel, 
                galaxy_density, 
                satellite_fraction
            )
    
    # Initial guess is the minimum point in the grid:
    i, j   = np.unravel_index( np.argmin(scoreGrid), scoreGrid.shape )
    x0, y0 = x[i, j], y[i, j]

    # Minimizing the score function to get the optimum values: NOTE: by default, search is done in the  
    # interval Mmin in [ 1e+12, 1e+13 ] and M1 in [ 1e+13, 1e+15 ].
    res = minimize(
        _scoreFunction, 
        x0     = [ x0, y0 ], 
        bounds = [
            ( np.log(1e+12), np.log(1e+13) ), # Mmin range
            ( np.log(1e+13), np.log(1e+15) ), # M1 range
        ], 
        args   = ( 
            halomodel, 
            galaxy_density, 
            satellite_fraction, 
         )
    )
    if not res.success:
        import warnings
        warnings.warn( f"Optimization failed with messagee {res.message!r} after {res.nit} iterations" )

    # Set the optimum values to the halo model
    score = _scoreFunction(
        res.x, 
        halomodel, 
        galaxy_density, 
        satellite_fraction, 
        return_model = True, 
    )
    print(f"Final score: {score:.4g}")
    print(f" - Galaxy density = {halomodel.averageGalaxyDensity():.3e} Mpc^-3")
    print(f" - Satellite fraction = {halomodel.averageSatelliteFraction():.3g}") 
    return halomodel



halomodel = _buildFromAbacusSimulation("AbacusSummit_hugebase_c000_ph000", 3.)
print( halomodel )