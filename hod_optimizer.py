# 
# Calculate optimum HOD parameters to get galaxy density and satellite fraction matching observation.
#

__version__ = "0.1a"

import logging
import numpy as np, numpy.typing
from scipy.optimize import minimize
from astropy.cosmology import w0waCDM
from galaxy_catalog_generator.halomodel import HaloModel

# Look-up table to get `sigma_8` values used in the abacus summit cosmologies. This data is not in the
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

def buildFromAbacusSimulation(
        simname      : str, 
        redshift     : float,
        sigmaM       : float = 0.0, 
        alpha        : float = 1.0, 
        powerspectrum: str   = "eisenstein98_zb", 
        massfunction : str   = "tinker08",  
    ) -> HaloModel:
    r"""
    Create a partially initialized halo model object using the cosmology parameters from the simulation. 
    Here ``M0`` is taken the same as ``Mmin``.

    Parameters
    ----------
    simname : str
        Name of the abacus simulation. 

    redshift : float
        Redshift for the simulation snapshot. This will raise error if the value is not available. 

    sigmaM : float default=0.0
        Width of the central galaxy transition range. A value of 0 means the relation is a step function.

    alpha : float, default=1.0
        Index for the  power law satelliete count relation.

    powerspectrum : str, default=`'eisenstein98_zb'`
        Power spectrum model to use. Parameters are taken from the simulation header.

    massfunction : str, default=`'tinker08'`
        Halo mass-function model to use. Parameters are taken from the simulation header.

    Returns
    -------
    halomodel : HaloModel
        A partially initialised halo model object, without ``Mmin``, ``M0`` and ``M1`` parameters. 
    
    """

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
        M0        = -1,       # using M0 = Mmin
        M1        = -1,       # M1 is calculated later...
        alpha     = alpha, 
        scaleSHMF = 0.5,      # values of sub-halo mass-function parameters do not affect other ...
        slopeSHMF = 2.0,      # HOD parameters: so a convenient value is set. 
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

def scoreFunction(
        x: tuple[float, float], 
        halomodel         : HaloModel,
        galaxy_density    : float, 
        satellite_fraction: float,
    ) -> float:
    r"""
    A estimate for the weighted difference of the calculated and required galaxy density and satellite 
    fraction values.

    Parameters
    ----------  
    x : tuple[float, float]

    halomodel : HaloModel
        Halo model object to use for calculations. All HOD parameters must be set.

    galaxy_density : float 
        Reference value for galaxy density at current redshift, in Mpc^-3 units.

    satellite_fraction : float
        Reference value for satellite fraction at current redshift.

    Returns
    -------
    score : float
        Log of the weighted difference. 

    """

    # Set the HOD paramters t the model
    Mmin, M1 = np.exp(x)
    halomodel._updateHaloParameters(Mmin = Mmin, M1 = M1, M0 = Mmin)

    # Relative diffrence in galaxy density
    deltaG = ( halomodel.averageGalaxyDensity( lnmb = np.log(1e+18) ) / galaxy_density - 1. )**2
    
    # Relative difference in satellite fraction
    deltaS = ( halomodel.averageSatelliteFraction( lnmb = np.log(1e+18) ) / satellite_fraction - 1. )**2
    
    # Total score: weighted distance from observed and calculated values
    score = np.log( deltaG + deltaS + 1e-16 ) 
    return score

def optimizeHaloModel(
        halomodel         : HaloModel,
        galaxy_density    : float, 
        satellite_fraction: float,
        Mmin_range        : tuple[float, float] = ( 1e+11, 1e+14 ),
        M1_range          : tuple[float, float] = ( 1e+13, 1e+15 ),
        gridsize          : int  = 12, 
    ) -> tuple[HaloModel, dict[str, numpy.typing.NDArray[np.float64]]]:
    r"""
    Optimize the HOD parameters ``Mmin`` and ``M1`` to match the calculated values of galaxy density and 
    satellite fraction with their required values.

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

    Returns
    -------
    halomodel : HaloModel
        Fully initialised halo model object with optimized HOD parameters. ``M0`` is same as ``Mmin``.

    datadict : dict
        A dict containing the values of score 

    """

    logger = logging.getLogger()

    # Creating the grid: this is used to locate the initial guess for the minimum
    logger.debug("generating grid...")
    ( xa, xb ), ( ya, yb ) = Mmin_range, M1_range
    x, y = np.meshgrid(
        np.linspace( np.log( xa ), np.log( xb ), gridsize ), # Mmin values
        np.linspace( np.log( ya ), np.log( yb ), gridsize ), # M1 values
    )
    scoreGrid = np.zeros_like(x)
    for i in range( scoreGrid.shape[0] ):
        for j in range( scoreGrid.shape[1] ):
            scoreGrid[i, j] = scoreFunction(
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
    logger.debug(f"optimizing for Mmin and M1 with guess ({x0=:.3f}, {y0=:.3f})...")
    logger.debug(f"using Mmin range [{xa:.3e}, {xb:.3e}] and M1 range [{ya:.3e}, {yb:.3e}]...")
    res = minimize(
        scoreFunction, 
        x0     = [ x0, y0 ], 
        bounds = [
            ( np.log(xa), np.log(xb) ), # Mmin range
            ( np.log(ya), np.log(yb) ), # M1 range
        ], 
        args   = ( 
            halomodel, 
            galaxy_density, 
            satellite_fraction, 
         )
    )
    if not res.success:
        logger.warning( f"optimization failed with messagee {res.message!r} after {res.nit} iterations" )

    # Set the optimum values to the halo model
    score = scoreFunction(
        res.x, 
        halomodel, 
        galaxy_density, 
        satellite_fraction, 
    )
    galaxyDensity     = halomodel.averageGalaxyDensity()
    satelliteFraction = halomodel.averageSatelliteFraction()
    logger.debug(f"optimum values: Mmin={halomodel.Mmin:.3e} Msun, M1={halomodel.M1:.3e} Msun")
    logger.debug(f"Final score: {score:.4g}")
    logger.debug(f"Galaxy density: {galaxyDensity:.3e} Mpc^-3")
    logger.debug(f"Satellite fraction: {satelliteFraction:.3g}") 

    datadict = {
        "galaxyDensity"    : galaxyDensity,
        "satelliteFraction": satelliteFraction,
        "scoreFunction"    : score,
        "scoreGrid"        : scoreGrid,
        "Mmin_Range"       : Mmin_range,
        "M1_Range"         : M1_range,
        "gridsize"         : gridsize,
    }
    return halomodel, datadict


if __name__ == "__main__":

    import click

    # Format of the output printed to stdout:
    _outputString = (
        "# -- Used input model parameters --\n"
        "# H0           : {H0:.8g}          \n"
        "# Tcmb0        : {Tcmb0:.8g}       \n"
        "# Om0          : {Om0:.8g}         \n"
        "# Ob0          : {Ob0:.8g}         \n"
        "# Ode0         : {Ode0:.8g}        \n"
        "# w0           : {w0:.8g}          \n"
        "# wa           : {wa:.8g}          \n"
        "# ns           : {ns:.8g}          \n"
        "# sigma8       : {sigma8:.8g}      \n"
        "# Delta        : {Delta:d}         \n"
        "# redshift     : {redshift:.8g}    \n"
        "# powerspectrum: {powerspectrum:s} \n" 
        "# massfunction : {massfunction:s}  \n"
        "# \n"
        "# -- Optimum HOD parameters: **copy-paste these values to the parameters file** -- \n"
        "Mmin  : {Mmin:.8g}   \n"  
        "M1    : {M1:.8g}     \n"
        "M0    : {M0:.8g}     \n"
        "sigmaM: {sigmaM:.8g} \n"
        "alpha : {alpha:.8g}  \n"
    )

    @click.command
    @click.version_option(__version__, message = "HOD Optimizer %(version)s") # Add --version
    @click.option("--simname"  , type = str  ,                              help = "Name of abacus simulation"      )
    @click.option("--redshift" , type = float,                              help = "Redshift of the snapshots"      )
    @click.option("--galdens"  , type = float,                              help = "target galaxy density in Mpc^-3")
    @click.option("--satfrac"  , type = float,                              help = "Target satellite fraction"      )
    @click.option("--sigma"    , type = float, default = 0.0              , help = "Central count transition width" )
    @click.option("--alpha"    , type = float, default = 1.0              , help = "Satellite count function slope" )
    @click.option("--powerspec", type = str  , default = "eisenstein98_zb", help = "Power spectrum model to use"    ) 
    @click.option("--massfunc" , type = str  , default = "tinker08"       , help = "Halo mass-function model to use")
    def _cli(
            simname  : str, 
            redshift : float,
            galdens  : float,
            satfrac  : float,
            sigma    : float = 0.0, 
            alpha    : float = 1.0, 
            powerspec: str   = "eisenstein98_zb", 
            massfunc : str   = "tinker08",      
        ) -> None:
        r"""
        Optimize the halo model parameters to get required galaxy density and satellite fraction. Copy paste 
        these values (printed to ``stdout``) to the parameter file for galaxy generator code to use this 
        model. 
        """

        assert  simname , "'simname'  is required"
        assert redshift, "'redshift' is required"
        assert galdens , "'galdens'  is required"
        assert satfrac , "'satfrac'  is required"

        halomodel, tree = optimizeHaloModel(
            halomodel = buildFromAbacusSimulation(
                simname       = simname, 
                redshift      = redshift, 
                sigmaM        = sigma, 
                alpha         = alpha, 
                powerspectrum = powerspec,
                massfunction  = massfunc,
            ),
            galaxy_density     = galdens, 
            satellite_fraction = satfrac, 
        )

        params = {
            "Mmin"         : halomodel.Mmin,
            "M1"           : halomodel.M1,
            "M0"           : halomodel.M0,
            "sigmaM"       : halomodel.sigmaM,
            "alpha"        : halomodel.alpha,
            "H0"           : halomodel.psmodel.cosmo.H0.value,
            "Tcmb0"        : halomodel.psmodel.cosmo.Tcmb0.value,
            "Om0"          : halomodel.psmodel.cosmo.Om0,
            "Ob0"          : halomodel.psmodel.cosmo.Ob0,
            "Ode0"         : halomodel.psmodel.cosmo.Ode0,
            "w0"           : halomodel.psmodel.cosmo.w0,
            "wa"           : halomodel.psmodel.cosmo.wa,
            "ns"           : halomodel.psmodel.ns,
            "sigma8"       : halomodel.psmodel.sigma8,
            "Delta"        : halomodel.Delta,
            "redshift"     : redshift,
            "powerspectrum": powerspec,
            "massfunction" : massfunc,
        }
        tree["optParameters"] = params

        # Save the data tree for later use
        import asdf
        af = asdf.AsdfFile(tree)
        af.write_to("hod_optimizer_result.asdf")
        af.close()

        # Write the parameter values to the stdout in a YML compatible format. All values except the HOD 
        # parameters are commented (this is for informative purpose only!) 
        print( _outputString.format( **params ) )

        return

    _cli()
