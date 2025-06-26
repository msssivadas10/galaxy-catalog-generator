#!/usr/bin/python3
# 
# CLI tools to process abacus halo catalogs and generate data.
# 

__version__ = "0.2a"

import os, os.path, glob, re, click, logging, asdf, yaml, numpy as np
from typing import IO, Literal, Any
from numpy.typing import NDArray
from astropy.cosmology import w0waCDM
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from galaxy_catalog_generator.halomodel import HaloModel

# Look-up table to get `sigma_8` values used in the abacus summit cosmologies. This data is not in the
# catalog headers, so it is taken from <https://abacussummit.readthedocs.io/en/latest/cosmologies.html>.
sigma8Table = {
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

############################################################################################################
#                                   HELPER FUNCTIONS FOR CALCULATIONS                                      #
############################################################################################################

def listCatalogs(simname: str, redshift: float, search_paths: list[str] = []) -> list[str]:
    r"""
    Search for abacus halo catalog files for given simulation name and redshift in the paths
    specified and return a list of paths to the files.  
    """

    logger = logging.getLogger()

    # Going through all search paths, looking for catalog files, until files are found...
    # If no search path specified, look in the current directory.
    globResult = []
    tailPath   = os.path.join(simname, "halos", f"z{redshift:.3f}", "halo_info", "halo_info_*.asdf")
    for path in search_paths:
        globResult = glob.glob( os.path.join(path, tailPath) )
        if globResult: break

    # Filter out only the files with names matchng the pattern and sorting based on the
    # integer index: 
    _files = { 
        int( m.group(1) ): m.string 
            for m in map( 
                lambda fn: re.search(r"halo_info_(\d{3}).asdf", fn), # for files like `halo_info_123.asdf`
                globResult 
            ) 
            if m 
    }
    files  = [ _files[idx] for idx in sorted(_files) ]

    if not files: 
        return logger.info(f"no files for simulation {simname!r} at redshift {redshift}: exiting...")
    logger.info(f"found {len(files)} files for simulation {simname!r} at redshift {redshift}")

    return files

def loadCatalog(fn: str) -> tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64]]:
    r"""
    Load data from a abacus halo catalog file. Halo ID, position in Mpc and mass in Msun are returned.
    """

    logger = logging.getLogger()
    
    logger.info(f"loading halo catalog from file: {fn!r}")
    catalog = CompaSOHaloCatalog(fn, cleaned = False, fields = ["id", "SO_central_particle", "N"])
    
    # Halo position coordinates in Mpc units
    hubble       = 0.01 * catalog.header["H0"]
    haloPosition = np.array( catalog.halos["SO_central_particle"] ) / hubble

    # Halo mass in Msun
    unitMass = catalog.header["ParticleMassMsun"]
    haloMass = np.array( catalog.halos["N"] ) * unitMass
    
    # Halo integer id
    haloID = np.array( catalog.halos["id"] )
    
    return haloID, haloPosition, haloMass

def loadCatalog1(files: list[str], Nrange: tuple[float, float] = (0., np.inf)) -> NDArray[np.float64]:
    r"""
    Load data from a abacus halo catalog file. Halo position in Mpc/h and mass in Msun/h 
    are returned. 
    
    This is different from ``loadCatalog``, since it loads data from all catalog files in a
    simulation, with an option for filtering based on halo mass. If no mass range is provided, 
    or the range is very large, the loaded data would also be large and may cause memory 
    issues.  
    """

    logger = logging.getLogger()
    
    Nmin, Nmax = Nrange
    retval     = []
    for fn in files:
        logger.info(f"loading halo catalog from file: {fn!r}, mass range: [{Nmin}, {Nmax}]")
        catalog = CompaSOHaloCatalog(fn, cleaned = False, fields  = ["SO_central_particle", "N"])
        
        # Halo mass in Msun/h
        unitMass = catalog.header["ParticleMassHMsun"]
        haloMass = np.array( catalog.halos["N"] )

        # -- Filtering based on mass
        mask     = ( haloMass >= Nmin ) & ( haloMass <= Nmax )
        haloMass = haloMass[mask] * unitMass
        
        # Halo position coordinates in Mpc/h units
        haloPosition = np.array( catalog.halos["SO_central_particle"] )[mask, :]

        retval.append( np.hstack([haloPosition, haloMass.reshape((-1, 1))]) )

    return np.vstack(retval)

def buildHaloModel(simname: str, redshift: float, **haloargs) -> HaloModel:
    r"""
    Initialize a halo model object using the cosmology parameters corresponding to the given abacus 
    simulation details and halo model parameters.
    """

    logger = logging.getLogger()

    import abacusnbody.metadata

    logger.info(f"loading metatdata for simulation {simname!r}for redshift {redshift}...")
    tree = abacusnbody.metadata.get_meta(simname, redshift)
    
    logger.info(f"creating halo model object...")
    haloargs = {
        "Mmin"     : None,
        "sigmaM"   : 0.  ,
        "M1"       : None,
        "M0"       : None,
        "alpha"    : 1.  ,
        "scaleSHMF": 0.5 ,
        "slopeSHMF": 2.  ,
        **haloargs
    }
    for __key in [ "redshift", "cosmo", "ns", "sigma8", "Delta", "meta" ]: 
        # These values are loaded from metadata...
        _ = haloargs.pop( __key, None )
    
    model = HaloModel.create(
        **haloargs,
        redshift  = tree["Redshift"], 
        cosmo     = w0waCDM(
            H0    = tree["H0"      ], 
            Om0   = tree["Omega_M" ], 
            Ode0  = tree["Omega_DE"],
            Ob0   = tree["omega_b" ] / ( 0.01*tree["H0"] )**2,
            w0    = tree["w0"      ], 
            wa    = tree["wa"      ], 
            Tcmb0 = 2.728, # CMB temperature value if fixed as 2.728 K
            name  = tree["SimName" ],
        ), 
        ns        = tree["n_s"], 
        sigma8    = sigma8Table[ re.search(r"c(\d{3})", tree["SimName"]).group(1) ], 
        Delta     = tree["SODensity"][0],
        meta      = {},
    )

    # Save a few other values from the simulation metadata. Only the values of ``simname`` and 
    # ``redshift`` are enough, beacuse using them one can get the full metadata using the package
    # ``abacusnbody.metadata``. 
    keymapping = { 
        "boxsize"     : "BoxSize"   , "particleMass"    : "ParticleMassHMsun", 
        "boxsizeMpc"  : "BoxSizeMpc", "particleMassMsun": "ParticleMassMsun" , 
    }
    model.meta.update({ altkey: tree[key] for altkey, key in keymapping.items() })  
    return model

def halomodelDict(model: HaloModel, include_meta: bool = True) -> dict[str, Any]:
    r"""
    Convert a halo model object into a mapping, ready to be serialized as YAML/JSON. 
    """

    cosmo = model.cosmo
    tree  = {
        "Mmin"     : model.Mmin, 
        "M0"       : model.M0, 
        "M1"       : model.M1, 
        "sigmaM"   : model.sigmaM, 
        "alpha"    : model.alpha,
        "scaleSHMF": model.scaleSHMF, 
        "slopeSHMF": model.slopeSHMF, 
        "Delta"    : model.Delta, 
        "ns"       : model.psmodel.ns,
        "sigma8"   : model.psmodel.sigma8,
        "H0"       : np.array(cosmo.H0).tolist(), 
        "Om0"      : cosmo.Om0, 
        "Ode0"     : cosmo.Ode0, 
        "Ob0"      : cosmo.Ob0, 
        "w0"       : cosmo.w0, 
        "wa"       : cosmo.wa, 
        "Tcmb0"    : np.array(cosmo.Tcmb0).tolist(),
        "simname"  : cosmo.name,
        "redshift" : model.redshift,
        "psmodel"  : model.psmodel.id,
        "mfmodel"  : model.mfmodel.id
    }
    if include_meta: 
        tree.update(model.meta)

    return tree

############################################################################################################
#                                       ARGUMENT VALIDATOR FUNCTIONS                                       #
############################################################################################################

def siminfoValidator(ctx: Any, param: Any, value: tuple[str, float]) -> tuple[str, float]:  
    r"""
    Check if the value is a tuple of valid abacus simulation name and redshift value.
    """

    try:
        import abacusnbody.metadata

        # Try to get the metadata for this simulation. This will raise an error if the details are 
        # incorrect...
        simname, redshift = value
        abacusnbody.metadata.get_meta(simname, redshift)
    except ValueError as e:
        raise click.BadParameter(e)  
    
    return value

def rangeValidator(ctx: Any, param: Any, value: tuple[float, float]) -> tuple[float, float]:  
    r"""
    Check if the value is a valid range specifier for a positive real number.
    """
      
    try:
        left, right = value
        assert left > 0. and right > 0, "values must be positive"
        assert left < right, "lower limit must be less than upper limit"
    except Exception as e:
        raise click.BadParameter(e)  
    
    return value

############################################################################################################
#                                               CLI AND COMMANDS                                           #
############################################################################################################

@click.group
@click.version_option(__version__, message = "Abacus Halos %(version)s")
def cli() -> None:
    r"""
    Tools using abacus simulation halo catalogs to generate data.
    """
    
    import logging.config, warnings
    
    warnings.catch_warnings(action="ignore")

    # Configure the loggers
    logPath = "logs"
    logFile = os.path.join( logPath, f"{ os.path.basename(__file__).rsplit('.', 1)[0] }.log" )
    os.makedirs( logPath, exist_ok = True )
    logging.config.dictConfig({
        "version": 1, 
        "disable_existing_loggers": True, 
        "formatters": {
            "default": { "format": "[ %(asctime)s %(levelname)s %(process)d ] %(message)s" }
        }, 
        "handlers": {
            "stream": {
                "level": "INFO", 
                "formatter": "default", 
                "class": "logging.StreamHandler", 
                "stream": "ext://sys.stdout"
            }, 
            "file": {
                "level": "INFO", 
                "formatter": "default", 
                "class": "logging.handlers.RotatingFileHandler", 
                "filename": logFile, 
                "mode": "a", 
                "maxBytes": 2097152, # create a new file if size exceeds 2 MiB
                "backupCount": 4     # use maximum 4 files
            }
        }, 
        "loggers": {
            "root": {
                "level": "INFO", 
                "handlers": [
                    "stream", 
                    "file"
                ]
            }
        }
    })
    return

@cli.command
@click.option("--siminfo", 
              type     = (str, float), 
              required = True, 
              callback = siminfoValidator, 
              help     = "A tuple of a valid abacus simulation name and redshift value", )
@click.option("-o", "--output-file", 
              type     = click.File("w"), 
              required = True, 
              help     = "Filename for the output (text CSV format)", )
@click.option("-m", "--mass-range", 
              type     = (float, float), 
              default  = ( 10., 1e+04 ), 
              callback = rangeValidator, 
              help     = "Values of the left- and right-most mass bins in particle mass unit", )
@click.option("--nbins", 
              type     = click.IntRange(min = 2), 
              default  = 32, 
              help     = "Number of bins", )
@click.option("-p", "--path" , 
              type     = click.Path(exists = True), 
              default  = [ os.getcwd() ], 
              multiple = True,  
              help     = "Paths to look for simulation catalog files", )
def massfunc(
        siminfo     : tuple[str, float], 
        output_file : IO, 
        mass_range  : tuple[float, float] = ( 10., 1e+04 ), 
        nbins       : int                 = 32, 
        path        : list[str]           = [], 
    ) -> None:
    r"""
    Estimate halo mass-function from catalog.

    Estimate halo mass-function from an abacus halo catalog. Estimated values will be saved as a CSV table
    of halo mass (Msun) and mass-function dn/dm (Mpc^-3), in natural log format.
    """

    logger = logging.getLogger()
    logger.info("command: massfunc")

    # List all catalog files in the given path:
    simname, redshift = siminfo
    files = listCatalogs(simname, redshift, path)
    if not files: return 

    # Loading particle mass and boxsize from the metadata:
    # NOTE: Halo model object is not needed for calculations, but created to get the simulation 
    # parameters, which will be written as the comments in the output file... 
    tree     = halomodelDict( model = buildHaloModel(simname, redshift) )
    unitMass = tree["particleMassMsun"] 
    boxsize  = tree["boxsizeMpc"]

    # Estimating halo mass-function
    massBinEdges = unitMass * np.logspace( np.log10(mass_range[0]), np.log10(mass_range[1]), nbins+1 ) 
    massFunc     = np.zeros(shape = massBinEdges.shape[0]-1)
    for file in files:
        _, _, haloMass = loadCatalog(file) # halo mass in Msun 
        massFunc[:]   += np.histogram(haloMass, bins = massBinEdges)[0][:]
    
    dlnMassH  = np.diff( np.log( massBinEdges ) )
    massH     = np.sqrt( massBinEdges[1:] * massBinEdges[:-1] )
    massFunc /= ( boxsize**3 * massH * dlnMassH ) # dn/dm in Mpc^-3

    nonzeroMask       = massFunc > 0.
    massFunctionTable = np.log( np.stack(( massH[nonzeroMask], massFunc[nonzeroMask] ), axis = -1) )
     
    # Save data as plain text CSV format
    logger.info(f"saving data to file: {output_file.name!r}")
    with output_file:

        # Header string (comment) contains simulation parameters for quick reference...
        tree = { 
            _key: tree[ _key ] for _key in [
                    "simname", "redshift", "particleMass", "boxsize", "H0", "Om0", "Ode0", "Ob0", 
                    "w0", "wa", "ns", "sigma8", "Tcmb0", "Delta",
            ] 
        }
        tree["table"] = {
            "columns": [
                { "name": "mass", "fmt": "log", "unit": "Msun"   },
                { "name": "dndm", "fmt": "log", "unit": "Mpc^-3" },
            ],
            "size": np.size(massFunctionTable, 0)
        }

        np.savetxt(
            output_file, 
            massFunctionTable, 
            header    = yaml.safe_dump(tree, explicit_start = True, explicit_end = True, sort_keys = False),
            delimiter = ',',
            comments  = '#'
        )
    return

@cli.command
@click.option("--siminfo", 
              type     = (str, float), 
              required = True, 
              callback = siminfoValidator, 
              help     = "A tuple of a valid abacus simulation name and redshift value", )
@click.option("--galaxy-density", "--ngal", 
              type     = click.FloatRange(min = 0., min_open = True), 
              required = True, 
              help     = "Observed value of the galaxy number density in Mpc^-3", )
@click.option("--satellite-fraction", "--fsat",  
              type     = click.FloatRange(0., 1.), 
              required = True, 
              help     = "Observed value of the satellite galaxy fraction", )
@click.option("-o", "--output-file", 
              type     = click.File("w"), 
              required = True, 
              help     = "Filename for the output (YAML format)", )
@click.option("--sigma-m", 
              type     = click.FloatRange(min = 0), 
              default  = 0., 
              help     = "Transition width for central galaxy count function", )
@click.option("--alpha", 
              type     = click.FloatRange(min = 0), 
              default  = 1., 
              help     = "Index for the power law satellite count function", )
@click.option("--scale-shmf", 
              type     = click.FloatRange(0., 1.), 
              default  = 0.5, 
              help     = "Scale for the SHMF - specify the maximum mass of a subhalo, as a fraction of halo", )
@click.option("--slope-shmf", 
              type     = click.FloatRange(min = 0., min_open = True), 
              default  = 2., 
              help     = "Slope for the SHMF", )
@click.option("--powerspectrum", "--mps", 
              type     = str, 
              default  = "eisenstein98_zb", 
              help     = "Matter power spectrum model used for calculations", )
@click.option("--massfunction", "--hmf", 
              type     = str, 
              default  = "tinker08", 
              help     = "Halo mass-function model used for calculations or path to the data file", )
@click.option("--mmin-range", 
              type     = (float, float), 
              default  = ( 1e+11, 1e+14 ), 
              callback = rangeValidator, 
              help     = "Search range for parameter Mmin (or M0) in Msun units", )
@click.option("--m1-range", 
              type     = (float, float), 
              default  = ( 1e+13, 1e+15 ), 
              callback = rangeValidator, 
              help     = "Search range for parameter M1 in Msun units", )
@click.option("--gridsize", 
              type     = click.IntRange(min = 4), 
              default  = 12, 
              help     = "Size of the Mmin-M1 grid used to guess minimum", )
def halomod(
        siminfo            : tuple[str, float], 
        galaxy_density     : float, 
        satellite_fraction : float,
        output_file        : IO, 
        sigma_m            : float                = 0. , 
        alpha              : float                = 1. ,
        scale_shmf         : float                = 0.5, 
        slope_shmf         : float                = 2. , 
        powerspectrum      : str                  = "eisenstein98_zb", 
        massfunction       : str                  = "tinker08", 
        mmin_range         : tuple[float, float]  = ( 1e+11, 1e+14 ),
        m1_range           : tuple[float, float]  = ( 1e+13, 1e+15 ),
        gridsize           : int                  = 12, 
    ) -> None:
    r"""
    Calculate halo model parameters.

    Calculate the optimum values for halo model, based on observed values of galaxy density and satellite 
    fraction. Other model parameters including cosmology are selected corresponding to the given abacus 
    simulation details.
    """

    from galaxy_catalog_generator.misc.halomodel_optimization import HODOptimizer

    logger = logging.getLogger()
    logger.info("command: halomod")

    # Creating the halo model object 
    simname, redshift = siminfo
    model = buildHaloModel(
        simname, 
        redshift, 
        Mmin      = None, 
        M1        = None, 
        M0        = None, 
        sigmaM    = sigma_m, 
        alpha     = alpha,
        scaleSHMF = scale_shmf, 
        slopeSHMF = slope_shmf, 
        psmodel   = powerspectrum, 
        mfmodel   = massfunction,
    )

    # Optimizing the model:
    logger.info(f"optimizing halo model...")
    optimizer = HODOptimizer(
        model, 
        galaxy_density     = galaxy_density, 
        satellite_fraction = satellite_fraction, 
        Mmin_range         = mmin_range, 
        M1_range           = m1_range, 
        gridsize           = gridsize, 
    )

    res = optimizer.optimize()
    logger.info(f"optimization exited at iteration {res.nit} with status={res.status} and message={res.message!r}")
    logger.info(f"optimization result: Mmin={res.Mmin:.4g} and M1={res.M1:.4g}")
    logger.info(f"optimum values: galaxy density={res.galaxy_density:.3g}, satellite fraction={res.satellite_fraction:.3g} (score={res.score:.3g})")

    # Saving file
    logger.info(f"saving model to {output_file.name!r}...")
    with output_file:
        yaml.safe_dump(
            {
                **halomodelDict(model, include_meta = 1), 
                "Settings" : {
                    "SimName"          : simname,
                    "Redshift"         : redshift,
                    "GalaxyDensity"    : galaxy_density, 
                    "SatelliteFraction": satellite_fraction,
                    "Mmin_Range"       : mmin_range, 
                    "M1_Range"         : m1_range
                },
                "OptimizationResult": {
                    "OptGalaxyDensity"    : float(res.galaxy_density), 
                    "OptSatelliteFraction": float(res.satellite_fraction),
                    "OptScore"            : float(res.score),
                    "Status"              : res.status, 
                    "Message"             : res.message,
                },
            }, 
            stream    = output_file, 
            sort_keys = False, 
        )
    return

@cli.command
@click.option("--siminfo", 
              type     = (str, float), 
              required = True, 
              callback = siminfoValidator, 
              help     = "A tuple of a valid abacus simulation name and redshift value", )
@click.option("--halomodel", 
              type     = click.File("r"), 
              required = True, 
              help     = "Path to the file which store the halo model parameters in YAML format", )
@click.option("--output-path", 
              type     = click.Path(file_okay = False), 
              required = True, 
              help     = "Path to the directory to save galaxy catalog files", )
@click.option("-p", "--path" , 
              type     = click.Path(exists = True), 
              default  = [ os.getcwd() ], 
              multiple = True,  
              help     = "Paths to look for simulation catalog files", )
def galaxies(
        siminfo     : tuple[str, float],
        halomodel   : IO, 
        output_path : str, 
        path        : list[str]  = [],
    ) -> None:
    r"""
    Galaxy catalog generation using halos.

    Generate a galaxy catalog using the specified abacus halo catalog and halo model parameters. 
    NOTE: Always make sure that the halo model parameters match the cosmology setup of the simulation 
    used and observation for more realistic results.
    """

    # Galaxy type ID: central galaxies are given the type number 1 and satellite galaxies 2 for  
    # easily identifying them from the catalog...
    CEN, SAT = 1, 2 

    logger = logging.getLogger()
    logger.info("command: galaxies")

    # List all catalog files in the given path:
    simname, redshift = siminfo
    files = listCatalogs(simname, redshift, path)
    if not files: return

    # Create the halo model object by loading the parameters from the given file
    haloargs = yaml.safe_load(halomodel)
    if not isinstance(haloargs, dict):
        return logger.error(f"cannot load halo model parameters from file {halomodel.name!r}")
    # This file can also contain cosmology model parameters. But, values taken from the simulation
    # metadata are used for consistency. NOTE: Make sure that the model parameters match the cosmology 
    # and redshift of the simulation. 
    haloargs  = { 
        key: haloargs[key] for key in haloargs
                            if key in  [
                                "Mmin", "M1", "M0", "sigmaM", "alpha", "scaleSHMF", "slopeSHMF", 
                                "psmodel", "mfmodel"
                            ]
    }

    try:
        model = buildHaloModel(simname, redshift, **haloargs)
    except Exception as e:
        return logger.error(f"error in creating halo model: {e}")
    
    # If the output path does not exist, create...
    output_path = os.path.join(output_path, simname, f"z{redshift:.3f}", "galaxy_info")
    os.makedirs(output_path, exist_ok = True)

    # Header: contains the values of cosmology and simulation parameters
    header = halomodelDict(model, include_meta = True)
        
    # Generate galaxy catalog for halos in each file
    filePattern = re.compile(r"halo_info_(\d+).asdf") # for files like `halo_info_123.asdf`
    for file in files:

        # Loading halo data
        # NOTE: position coordinates are in Mpc and mass are in Msun unit 
        haloID, haloPosition, haloMass = loadCatalog(file)
        
        # Placing central and satellite galaxies in each halo randomly, based on the halo model.
        # NOTE: This is done in loop - slow for large catalogs 
        logger.info(f"generating galaxy catalog for { len(haloID) } halos...")
        halosPerFile, cgalaxyPerFile, sgalaxyPerFile          = 0, 0, 0
        parentHaloID, galaxyPositions, galaxyMass, galaxyType = [], [], [], []
        for hid, posH, massH in zip(haloID, haloPosition, haloMass):
            halosPerFile += 1
            galaxyData    = model.generateSatellitePositions( np.log(massH), posH )
            if galaxyData.shape[0] < 1:
                continue
            cgalaxyPerFile += 1
            sgalaxyPerFile += galaxyData.shape[0]-1
            parentHaloID.append( np.repeat(hid, repeats = galaxyData.shape[0]) )
            galaxyPositions.append( galaxyData[:, 0:3] )
            galaxyMass.append( galaxyData[:, 3] )
            galaxyType.append(np.array( [CEN] + [SAT] * (galaxyData.shape[0]-1), dtype = np.uint8 ))

        galaxyPercentage = 100 * ( cgalaxyPerFile / halosPerFile )
        logger.info(f"placed galaxies in {cgalaxyPerFile} out of {halosPerFile} halos ({galaxyPercentage:.3f}%)")

        # NOTE: position and mass are saved in Mpc/h and Msun/h units, respectrively.
        parentHaloID    = np.hstack(parentHaloID)
        galaxyPositions = np.vstack(galaxyPositions) * model.cosmo.h # galaxy position in Mpc/h
        galaxyMass      = np.hstack(galaxyMass)      * model.cosmo.h # galaxy mass in Msun/h
        galaxyType      = np.hstack(galaxyType)
        
        # Wrapping around the coordinates that falls outside the simulation boxsize, with a 
        # periodic boundary condition, as in the original simulation...
        boxsize = model.meta["boxsize"] # in Mpc/h
        galaxyPositions = ( galaxyPositions + 0.5*boxsize ) % boxsize - 0.5*boxsize # in Mpc/h

        # Saving the catalog to a file in output path. File index is same as the halo file index
        outputFile = os.path.join(
            output_path, 
            f"galaxy_info_{ filePattern.match(os.path.basename(file)).group(1) }.asdf", 
        )
        logger.info(f"saving catalog to file: {outputFile!r}") 
        af = asdf.AsdfFile({
            "header"    : {
                **header, 
                "halosProcessed"   : halosPerFile, 
                "centralGalaxies"  : cgalaxyPerFile, 
                "satelliteGalaxies": sgalaxyPerFile,
            }, 
            "data"      : {
                "parentHaloID"  : parentHaloID, 
                "galaxyPosition": galaxyPositions, 
                "galaxyMass"    : galaxyMass, 
                "galaxyType"    : galaxyType,
            }
        })
        af.write_to(outputFile, all_array_compression = 'zlib')
        af.close()
        
    logger.info("galaxy catalog generation completed! :)")
    return

@cli.command
@click.option("--siminfo", 
              type     = (str, float), 
              required = True, 
              callback = siminfoValidator, 
              help     = "A tuple of a valid abacus simulation name and redshift value", )
@click.option("-o", "--output-file", 
              type     = click.File("wb"), 
              required = True, 
              help     = "Filename for the output (text CSV format)", )
@click.option("-m1", "--mrange1", 
              type     = (float, float), 
              default  = (0., np.inf),
              callback = rangeValidator, 
              help     = "Mass interval #1 in particle mass unit, for sub-catalog", )
@click.option("-m2", "--mrange2", 
              type     = (float, float), 
              default  = (0., np.inf),
              callback = rangeValidator, 
              help     = "Mass interval #2 in particle mass unit, for sub-catalog", )
@click.option("--estimator", 
              type       = click.Choice(["nat", "ls", "ham", "dp"], case_sensitive = False), 
              default    = "nat", 
              help       = "Specify the estimator to use", )
@click.option("-r", "--rrange", 
              type     = click.Tuple((float, float)), 
              required = True,
              callback = rangeValidator, 
              help     = "Distance bins range in Mpc/h", )
@click.option("--rbins", 
              type     = click.IntRange(min = 2), 
              default  = 32, 
              help     = "Number of distance bins", )
@click.option("--nthreads", 
              type     = click.IntRange(min = 0), 
              default  = 0, 
              help     = "Number of threads to use (0 for using all available threads)", )
@click.option("-p", "--path" , 
              type     = click.Path(exists = True), 
              default  = [ os.getcwd() ], 
              multiple = True,  
              help     = "Paths to look for simulation catalog files", )
def corrfunc(
        siminfo     : tuple[str, float], 
        output_file : IO, 
        mrange1     : tuple[float, float],
        mrange2     : tuple[float, float],
        rrange      : tuple[float, float],
        estimator  : Literal["nat", "ls", "ham", "dp"] = "nat",
        rbins       : int       =  16, 
        nthreads    : int       =  0,
        path        : list[str] =  [],
    ) -> None:
    r"""
    Estimate halo correlation function from catalog.

    Estimate halo 2-point correlation function (auto and cross) from an abacus halo catalog. 
    Estimated values will be saved as an asdf file.
    """

    from Corrfunc.theory.DD import DD
    # from Corrfunc.utils import convert_3d_counts_to_cf

    logger = logging.getLogger()
    logger.info("command: corrfunc")

    # List all catalog files in the given path:
    simname, redshift = siminfo
    files = listCatalogs(simname, redshift, path)
    if not files: return 

    # Loading catalogs: if the ranges are same, auto correlation is calculated and single
    # catalog is needed 
    autocorr = np.allclose(mrange1, mrange2)
    D1 = loadCatalog1(files, Nrange = mrange1)
    D2 = loadCatalog1(files, Nrange = mrange2) if not autocorr else D1

    # Generating random points (using R2 same as R1)
    import abacusnbody.metadata
    tree    = abacusnbody.metadata.get_meta(simname, redshift)
    boxsize = tree["BoxSize"] 

    Nr = 3*max( D1.shape[0], D2.shape[0] )

    logger.info(f"generating random catalog of size {Nr}...")
    R1 = np.random.uniform(0., boxsize, size = [Nr, 3])

    # Calculating correlation function
    rBinEdges = np.logspace( np.log10(rrange[0]), np.log10(rrange[1]), rbins+1 )
    rCenter   = np.sqrt( rBinEdges[1:] * rBinEdges[:-1] ) 

    # -- Calculating pair counts D1D2, D1R2, D2R1 and R1R2 for LS estimator:
    kwargs = dict(
        nthreads = nthreads or os.cpu_count(), 
        binfile  = rBinEdges, 
        periodic = True, 
        boxsize  = boxsize,
        verbose  = False,
    ) # common args to all pair count calls

    # Calculating number of pairs between D1 and R2 catalogs, D1R2, as it is common if D2 is
    # same as D1 or not...  
    logger.info(f"counting pairs of D1 and R...")
    D1R2 = DD(
        autocorr = 0, 
        X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
        X2 = R1[:,0], Y2 = R1[:,1], Z2 = R1[:,2], 
        **kwargs, 
    )
    D1R2 = D1R2["npairs"]
    
    if autocorr:
        # Here, D2 catalog is same as D1 catalog, so D1D2 is the number of pairs in the D1
        # catalog only (``autocorr = 1``) and D2R1 is same as D1R2... 
        logger.info(f"counting pairs of D1 and D1...")
        D1D2 = DD(
            autocorr = 1, 
            X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
            **kwargs, 
        )
        D1D2 = D1D2["npairs"]

        D2R1 = D1R2
    else:
        # D1 and D2 are different, calculating D1D2 and D2R1... 
        logger.info(f"counting pairs of D1 and D2...")
        D1D2 = DD(
            autocorr = 0, 
            X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
            X2 = D2[:,0], Y2 = D2[:,1], Z2 = D2[:,2], 
            **kwargs, 
        )     
        D1D2 = D1D2["npairs"]
    
        logger.info(f"counting pairs of D2 and R...")
        D2R1 = DD(
            autocorr = 0, 
            X1 = R1[:,0], Y1 = R1[:,1], Z1 = R1[:,2], 
            X2 = D2[:,0], Y2 = D2[:,1], Z2 = D2[:,2], 
            **kwargs, 
        )
        D2R1 = D2R1["npairs"]
    
    logger.info(f"counting pairs of R and R...")
    R1R2 = DD(
        autocorr = 1, 
        X1 = R1[:,0], Y1 = R1[:,1], Z1 = R1[:,2], 
        **kwargs, 
    )
    R1R2 = R1R2["npairs"]
    
    # -- Normalized pair counts:
    d1d2 = D1D2 / ( D1.shape[0]*D2.shape[0] )
    d1r2 = D1R2 / ( D1.shape[0]*R1.shape[0] )
    d2r1 = D2R1 / ( D2.shape[0]*R1.shape[0] )
    r1r2 = R1R2 / ( R1.shape[0]*R1.shape[0] )

    # -- Correlation function:
    estimator     = estimator.lower()
    estimatorName = {
        "nat" : "Natural", 
        "dp"  : "Davis-Peebles", 
        "ham" : "Hamilton", 
        "ls"  : "Landy-Szalay", 
    }.get(estimator, "Landy-Szalay")
    
    logger.info(f"calculating correlation function using {estimatorName!r} estimator...")
    if estimator == "nat":
        # Natural estimator
        xi = d1d2 / r1r2 - 1.
    elif estimator == "dp":
        # Davis-Peebles estimator
        xi = 2*d1d2 / (d1r2 + d2r1) + 1.
    elif estimator == "ham":
        # Hamilton estimator
        xi = (d1d2 / d1r2) * (r1r2 / d2r1) - 1.
    else:
        # Landy-Szalay estimator (default)
        xi = (d1d2 - d1r2 - d2r1 + r1r2) / r1r2
    # xi = convert_3d_counts_to_cf(D1.shape[0], D2.shape[0], Nr, Nr, D1D2, D1R2, D2R1, R1R2, estimator = 'LS')

    # Saving the output file:
    logger.info(f"saving correlation values to file {output_file.name!r}")
    af = asdf.AsdfFile({
        "header" : {
            k: tree[k] for k in [ 
                "SimName", "Redshift", "H0", "Omega_M", "Omega_DE", "omega_b", 
                "w0", "wa", "n_s", "SODensity", "BoxSize", "ParticleMassHMsun", 
                "BoxSizeMpc", "ParticleMassMsun",
            ]
        } | {
            "sigma8"       : sigma8Table[ re.search(r"c(\d{3})", tree["SimName"]).group(1) ], 
            "NRange1"      : list( mrange1 ), 
            "NRange2"      : list( mrange2 ), 
            "CatalogSize1" : D1.shape[0], 
            "CatalogSize2" : D2.shape[0],
            "RandomSize"   : R1.shape[0],
            "Estimator"    : estimatorName,  
        }, 
        "data" : {
            "r"    : rCenter, 
            "xi"   : xi,
            "D1D2" : D1D2, 
            "D1R2" : D1R2,
            "D2R1" : D2R1,
            "R1R2" : R1R2, 
        }
    })
    af.write_to(output_file, all_array_compression = "zlib")
    af.close()

    logger.info("correlation calculation completed!...")
    return

if __name__ == "__main__":
    cli()
    