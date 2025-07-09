#!/usr/bin/python3
# 
# CLI tools to process abacus halo catalogs and generate data.
# 

__version__ = "0.2b"

import os, os.path, glob, re, click, logging, asdf, yaml, numpy as np
from typing import IO, Literal, Any
from numpy.typing import NDArray
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from galaxy_catalog_generator.halomodel import HaloModel
from galaxy_catalog_generator.misc.correlation import PairCountData

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

catarray = np.dtype([("id", np.uint64), ("pos", np.float64, [3]), ("mass", np.float64)])

def loadCatalog(fn: str | list[str], Nrange: tuple[float, float] = None) -> NDArray[np.void]:
    r"""
    Load data from a abacus halo catalog file. Halo ID, position in Mpc and mass in Msun are returned.
    """

    logger = logging.getLogger()

    if isinstance(fn, (list, tuple)):
        # Load data from multiple files. (NOTE: loading a large data may cause memory issues...)
        return np.hstack([ loadCatalog(_fn, Nrange) for _fn in fn ])

    # Load data from a single file
    logger.info(f"loading halo catalog from file: {fn!r}")
    catalog = CompaSOHaloCatalog(fn, cleaned = False, fields = ["id", "SO_central_particle", "N"])

    # Halo integer id
    haloID = np.array( catalog.halos["id"] )
    data   = np.zeros(shape = len(haloID), dtype = catarray)
    data["id"] = haloID

    # Halo mass in Msun
    unitMass = catalog.header["ParticleMassMsun"]
    data["mass"] = np.array( catalog.halos["N"] ) * unitMass
    
    # Halo position coordinates in Mpc units
    hubble      = 0.01 * catalog.header["H0"]
    data["pos"] = np.array( catalog.halos["SO_central_particle"] ) / hubble
    
    # Filtering based on mass (optional)
    if Nrange is not None:
        logger.info(f"applying mass range {Nrange}...")
        data = data[
            ( data["mass"] >= Nrange[0] * unitMass ) & 
            ( data["mass"] <= Nrange[1] * unitMass )
        ]

    return data

def loadMetadata(simname: str, redshift: float, return_hm : bool = True, **kwargs) -> tuple[dict, HaloModel]:
    r"""
    Load a subset of the metadata corresponding to a simulation. Also initialise a halomodel object using 
    the loaded data. Other specified keywords arguments are passed to the object constructor. 
    """

    logger = logging.getLogger()
    logger.info(f"loading metatdata for simulation {simname!r}for redshift {redshift}...")

    import abacusnbody.metadata
    meta = abacusnbody.metadata.get_meta(simname, redshift)

    # Filter the dict to contain only the selected fields
    fields = [("SimName"          , "simname"          ), 
              ("Redshift"         , "redshift"         ), 
              ("H0"               , "H0"               ), 
              ("Omega_M"          , "Om0"              ), 
              ("Omega_DE"         , "Ode0"             ), 
              ("omega_b"          , "Ob0"              ), 
              ("w0"               , "w0"               ), 
              ("wa"               , "wa"               ), 
              ("n_s"              , "ns"               ), 
              ("SODensity"        , "SODensity"        ), 
              ("BoxSize"          , "boxsize"          ), 
              ("ParticleMassHMsun", "particleMassHMsun"), 
              ("BoxSizeMpc"       , "boxsizeMpc"       ), 
              ("ParticleMassMsun" , "particleMassMsun" ),]
    tree = { kn : meta[ko] for ko, kn in fields }

    tree["Ob0"   ] = tree["Ob0" ] / ( 0.01*tree["H0"] )**2 # converting Obh2 to Ob0
    tree["sigma8"] = sigma8Table[ re.search(r"c(\d{3})", simname).group(1) ] # sigma_8 parameter
    tree["Tcmb0" ] = 2.728 # CMB temperature value 

    if not return_hm:
        return tree, None

    # Creating the halo model object:
    from astropy.cosmology import w0waCDM

    # Matter power spectrum model
    psmodel = kwargs.get( "psmodel", "eisenstein98_zb" )
    tree["psmodel"] = psmodel

    # Halo mass-function model
    mfmodel = kwargs.get( "mfmodel", "tinker08" )
    tree["mfmodel"] = mfmodel

    mfFilePath = os.path.abspath(mfmodel)
    if os.path.isfile(mfFilePath):
        # Load mass function data from a file
        logger.info(f"found mass-function data file {mfFilePath!r}")
        
        if re.match(r".+\.asdf$", mfmodel): 
            # Data file is in ASDF format: mass-function data in a `data` section with fields 
            # `ln_M` (i.e., ln(M), with mass M in Msun units) and `ln_dndM` (i.e, ln(dndM), with
            # mass-function dndM in Mpc^-3 unit): 
            with asdf.open(mfmodel) as af:
                hmfdata = np.array([af["data"]["ln_M"   ], 
                                    af["data"]["ln_dndM"]]).T
        else:
            # Data file is a plain text file, comma seperated, with first two columns identified 
            # as the fields `ln_M` and `ln_dndM`: 
            hmfdata = np.loadtxt(mfmodel, usecols = [0, 1], delimiter = ',', comments = '#')

        from galaxy_catalog_generator.halomassfunction import massFunctionFromData
        
        logger.info("creating custom mass-function model from loaded data")
        mfmodel         = massFunctionFromData("data", hmfdata)
        mfmodel.id      = mfFilePath # change the model id to the file path 
        tree["mfmodel"] = mfFilePath

    logger.info(f"creating halo model object...")
    cargs = ( "H0", "Om0", "Ode0", "Ob0" ,"w0", "wa", "Tcmb0" ) 
    model = HaloModel.create(
        Mmin      = kwargs.get( "Mmin"     , -1.  ),
        M0        = kwargs.get( "M0"       , -1.  ),
        M1        = kwargs.get( "M1"       , -1.  ),
        sigmaM    = kwargs.get( "sigmaM"   ,  0.  ),
        alpha     = kwargs.get( "alpha"    ,  2.  ),
        scaleSHMF = kwargs.get( "scaleSHMF",  0.5 ),
        slopeSHMF = kwargs.get( "slopeSHMF",  2.  ),
        redshift  = tree["redshift"],
        cosmo     = w0waCDM( name = tree["simname" ], **{ k : tree[k] for k in cargs } ), 
        psmodel   = psmodel,
        mfmodel   = mfmodel,
        ns        = tree["ns"],
        sigma8    = tree["sigma8"],
        Delta     = tree["SODensity"][0],
        meta      = None
    )

    attrs = [ "Mmin", "sigmaM", "M0", "M1", "alpha", "scaleSHMF", "slopeSHMF", ] # HOD parameters
    tree.update({ attr: getattr(model, attr) for attr in attrs })

    return tree, model

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
    tree, _  = loadMetadata(simname, redshift, return_hm = False)
    tree.update({ "Nrange" : mass_range, "nbins" : nbins })
    unitMass = tree["particleMassMsun"] 
    boxsize  = tree["boxsizeMpc"]

    # Estimating halo mass-function
    massBinEdges = unitMass * np.logspace( np.log10(mass_range[0]), np.log10(mass_range[1]), nbins+1 ) 
    haloCount    = np.zeros(shape = massBinEdges.shape[0]-1)
    for file in files:
        data = loadCatalog(file) # halo mass in Msun 
        haloCount[:]  += np.histogram(data["mass"], bins = massBinEdges)[0][:]
    
    dlnMassH = np.log( massBinEdges[1] ) - np.log( massBinEdges[0] ) 
    massH    = np.sqrt( massBinEdges[1:] * massBinEdges[:-1] )

    # Filtering non zero values
    nonzeroMask  = np.argwhere(haloCount > 0.).flatten()
    massH, haloCount = massH[nonzeroMask], haloCount[nonzeroMask]

    massFunc    = haloCount / ( boxsize**3 * massH * dlnMassH ) # dn/dm in Mpc^-3
    massFuncErr = np.sqrt(haloCount) / haloCount # relative error assuming poisson noise
     
    # Save data as asdf file
    logger.info(f"saving data to file: {output_file.name!r}")
    with asdf.AsdfFile({
            "header" : tree, 
            "data"   : {
                "ln_M"     : np.log(massH),
                "ln_dndM"  : np.log(massFunc),
                "relError" : massFuncErr,
            }
        }) as af:
        af.write_to(output_file)

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
    tree, model = loadMetadata(
        simname, 
        redshift, 
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
    tree.update({ "Mmin": model.Mmin, "M0": model.M0, "M1": model.M1 })

    # Saving file
    logger.info(f"saving model to {output_file.name!r}...")
    with output_file:
        yaml.safe_dump(
            {
                **tree, 
                "settings" : {
                    "galaxyDensity"     : galaxy_density, 
                    "satelliteFraction" : satellite_fraction,
                    "Mmin_Range"        : mmin_range, 
                    "M1_Range"          : m1_range
                },
                "optimizationResult": {
                    "optGalaxyDensity"     : float(res.galaxy_density), 
                    "optSatelliteFraction" : float(res.satellite_fraction),
                    "optScore"             : float(res.score),
                    "status"               : res.status, 
                    "message"              : res.message,
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
    with halomodel:
        haloargs = yaml.safe_load(halomodel)
    if not isinstance(haloargs, dict):
        return logger.error(f"cannot load halo model parameters from file {halomodel.name!r}")
    # This file can also contain cosmology model parameters. But, values taken from the simulation
    # metadata are used for consistency. NOTE: Make sure that the model parameters match the cosmology 
    # and redshift of the simulation. 
    haloargs  = { 
        key: haloargs[key] for key in [
            "Mmin", "M1", "M0", "sigmaM", "alpha", "scaleSHMF", "slopeSHMF", 
            "psmodel", "mfmodel"
        ]
    }

    try:
        tree, model = loadMetadata(simname, redshift, **haloargs)
    except Exception as e:
        return logger.error(f"error in creating halo model: {e}")
    
    # If the output path does not exist, create...
    output_path = os.path.join(output_path, simname, f"z{redshift:.3f}", "galaxy_info")
    os.makedirs(output_path, exist_ok = True)

    # Generate galaxy catalog for halos in each file
    filePattern = re.compile(r"halo_info_(\d+).asdf") # for files like `halo_info_123.asdf`
    for file in files:

        # Loading halo data
        # NOTE: position coordinates are in Mpc and mass are in Msun unit 
        data = loadCatalog(file)
        haloID, haloPosition, haloMass = data["id"], data["pos"], data["mass"]
        
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
        # TODO: save in Mpc and Msun
        parentHaloID    = np.hstack(parentHaloID)
        galaxyPositions = np.vstack(galaxyPositions) * model.cosmo.h # galaxy position in Mpc/h
        galaxyMass      = np.hstack(galaxyMass)      * model.cosmo.h # galaxy mass in Msun/h
        galaxyType      = np.hstack(galaxyType)
        
        # Wrapping around the coordinates that falls outside the simulation boxsize, with a 
        # periodic boundary condition, as in the original simulation...
        boxsize = tree["boxsize"] # in Mpc/h
        galaxyPositions = ( galaxyPositions + 0.5*boxsize ) % boxsize - 0.5*boxsize # in Mpc/h

        # Saving the catalog to a file in output path. File index is same as the halo file index
        outputFile = os.path.join(
            output_path, 
            f"galaxy_info_{ filePattern.match(os.path.basename(file)).group(1) }.asdf", 
        )
        logger.info(f"saving catalog to file: {outputFile!r}") 
        with asdf.AsdfFile({
                "header" : {
                    **tree, 
                    "halosProcessed"    : halosPerFile, 
                    "centralGalaxies"   : cgalaxyPerFile, 
                    "satelliteGalaxies" : sgalaxyPerFile,
                }, 
                "data" : {
                    "parentHaloID"   : parentHaloID, 
                    "galaxyPosition" : galaxyPositions, 
                    "galaxyMass"     : galaxyMass, 
                    "galaxyType"     : galaxyType,
                }
            }) as af:
            af.write_to(outputFile, all_array_compression = 'zlib')
        
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
        estimator   : Literal["nat", "ls", "ham", "dp"] = "nat",
        rbins       : int       =  16, 
        nthreads    : int       =  0,
        path        : list[str] =  [],
    ) -> None:
    r"""
    Estimate halo correlation function from catalog.

    Estimate halo 2-point correlation function (auto and cross) from an abacus halo catalog. 
    Estimated values will be saved as an asdf file.
    """

    logger = logging.getLogger()
    logger.info("command: corrfunc")

    # List all catalog files in the given path:
    simname, redshift = siminfo
    files = listCatalogs(simname, redshift, path)
    if not files: return 
    
    tree, _ = loadMetadata(simname, redshift, return_hm = False)
    hubble  = tree["H0"] / 100.

    # Loading catalogs: if the ranges are same, auto correlation is calculated and single
    # catalog is needed 
    autocorr = np.allclose(mrange1, mrange2)
    D1 = loadCatalog(files, Nrange = mrange1)["pos"]
    D2 = loadCatalog(files, Nrange = mrange2)["pos"] if not autocorr else None

    # Calculating correlation function
    estimator = estimator.lower()
    rBinEdges = np.logspace( np.log10(rrange[0]), np.log10(rrange[1]), rbins+1 ) # in Mpc/h
    rCenter   = np.sqrt( rBinEdges[1:] * rBinEdges[:-1] ) 
    countData = PairCountData.countPairs(
        D1, 
        D2, 
        rBinEdges / hubble, # in Mpc
        boxsize  = tree["boxsizeMpc"], 
        periodic = True, 
        diff_r   = False, 
        nthreads = nthreads or os.cpu_count(),
    )
    xi, xiError = countData.correlation(estimator)
    tree.update({
        "NRange1"   : list( mrange1 ), 
        "NRange2"   : list( mrange2 ), 
        "sizeD1"    : countData.ND1, 
        "sizeD2"    : countData.ND2,
        "sizeR"     : countData.NR1,
        "estimator" : estimator,  
    })

    # Saving the output file:
    logger.info(f"saving correlation values to file {output_file.name!r}")
    with asdf.AsdfFile({
            "header" : tree, 
            "data"   : {
                "r"     : rCenter, 
                "xi"    : xi,
                "xiErr" : xiError,
                "D1D2"  : countData.D1D2, 
                "D1R2"  : countData.D1R2,
                "D2R1"  : countData.D2R1,
                "R1R2"  : countData.R1R2, 
            }
        }) as af:
        af.write_to(output_file, all_array_compression = "zlib")

    logger.info("correlation calculation completed!...")
    return

if __name__ == "__main__":
    cli()
    