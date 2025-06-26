#!/usr/bin/python3
# 
# CLI tools to calculate cross correlations from abacus halo catalogs.
# 

__version__ = "0.1a"

import os, os.path, glob, re, click, logging, numpy as np
from typing import Any, IO
from numpy.typing import NDArray
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from Corrfunc.theory.DD import DD
from Corrfunc.utils import convert_3d_counts_to_cf

# Look-up table to get `sigma_8` values used in the abacus summit cosmologies. This data is not in the
# catalog headers, so it is taken from <https://abacussummit.readthedocs.io/en/latest/cosmologies.html>.
AbacusCosmologySigma8Table = {
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

def loadCatalog(files: list[str], Nmin: float = 0., Nmax: float = np.inf) -> NDArray[np.float64]:
    r"""
    Load data from a abacus halo catalog file. Halo position in Mpc/h and mass in Msun/h 
    are returned. The catalog can be optionally filtered to only include halos with mass in 
    a specified range
    """

    logger = logging.getLogger()
    
    retval = []
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
    Check if the value is a valid distance range specifier.
    """

    try:
        left, right = value
        assert left > 0. and right > 0, "distance values must be positive"
        assert left < right, "lower limit must be less than upper limit"
    except Exception as e:
        raise click.BadParameter(e)  
    
    return value

############################################################################################################
#                                               CLI AND COMMANDS                                           #
############################################################################################################

@click.command
@click.version_option(__version__, message = "Abacus Correlation %(version)s")
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
def corrcalc(
        siminfo     : tuple[str, float], 
        output_file : IO, 
        mrange1     : tuple[float, float],
        mrange2     : tuple[float, float],
        rrange      : tuple[float, float],
        rbins       : int       = 16, 
        nthreads    : int       = 0,
        path        : list[str] = [],
    ) -> None:
    r"""
    Estimate halo correlation function from catalog.

    Estimate halo 2-point correlation function (auto and cross) from an abacus halo catalog. 
    Estimated values will be saved as an asdf file.
    """

    import warnings
    warnings.catch_warnings(action="ignore")
    
    # Configure the loggers
    import logging.config

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

    # ==================================================================================================== #
    
    logger = logging.getLogger()

    # List all catalog files in the given path:
    simname, redshift = siminfo
    files = listCatalogs(simname, redshift, path)
    if not files: return 

    # Loading catalogs: if the ranges are same, auto correlation is calculated and single
    # catalog is needed 
    autocorr = np.allclose(mrange1, mrange2)
    D1 = loadCatalog(files, Nmin = mrange1[0], Nmax = mrange1[1])
    D2 = loadCatalog(files, Nmin = mrange2[0], Nmax = mrange2[1]) if not autocorr else D1

    # Generating random points
    import abacusnbody.metadata
    tree    = abacusnbody.metadata.get_meta(simname, redshift)
    boxsize = tree["BoxSize"] 

    Nr = 3*max( D1.shape[0], D2.shape[0] )

    logger.info(f"generating random catalog of size {Nr}...")
    R1 = np.random.uniform(0., boxsize, size = [Nr, 3])

    # Calculating correlation function
    rBinEdges = np.linspace(rrange[0], rrange[1], rbins+1)
    rCenter   = np.sqrt( rBinEdges[1:] * rBinEdges[:-1] ) 

    # -- Calculating pair counts D1D2, D1R, D2R and RR for LS estimator:
    kwargs = dict(
        nthreads = nthreads or os.cpu_count(), 
        binfile  = rBinEdges, 
        periodic = True, 
        boxsize  = boxsize,
        verbose  = False,
    ) # common args to all pair count calls

    # Calculating number of pairs between D1 and R catalogs, D1R, as it is common if D2 is
    # same as D1 or not...  
    logger.info(f"counting pairs of D1 and R...")
    D1R  = DD(
        autocorr = 0, 
        X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
        X2 = R1[:,0], Y2 = R1[:,1], Z2 = R1[:,2], 
        **kwargs, 
    )
    
    if autocorr:
        # Here, D2 catalog is same as D1 catalog, so D1D2 is the number of pairs in the D1
        # catalog only (``autocorr = 1``) and D2R is same as D1R... 
        logger.info(f"counting pairs of D1 and D1...")
        D1D2 = DD(
            autocorr = 1, 
            X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
            **kwargs, 
        )

        D2R  = D1R
    else:
        # D1 and D2 are different, calculating D1D2 and D2R... 
        logger.info(f"counting pairs of D1 and D2...")
        D1D2 = DD(
            autocorr = 0, 
            X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
            X2 = D2[:,0], Y2 = D2[:,1], Z2 = D2[:,2], 
            **kwargs, 
        )     
    
        logger.info(f"counting pairs of D2 and R...")
        D2R  = DD(
            autocorr = 0, 
            X1 = R1[:,0], Y1 = R1[:,1], Z1 = R1[:,2], 
            X2 = D2[:,0], Y2 = D2[:,1], Z2 = D2[:,2], 
            **kwargs, 
        )
    
    logger.info(f"counting pairs of R and R...")
    RR   = DD(
        autocorr = 1, 
        X1 = R1[:,0], Y1 = R1[:,1], Z1 = R1[:,2], 
        **kwargs, 
    )
    
    # -- Correlation function:
    logger.info(f"calculating correlation function (LS estimator)...")
    xi = convert_3d_counts_to_cf(D1.shape[0], D2.shape[0], Nr, Nr, D1D2, D1R, D2R, RR, estimator = 'LS')

    # Saving the output file:
    import asdf

    logger.info(f"saving correlation values to file {output_file.name!r}")
    af = asdf.AsdfFile({
        "header" : {
            k: tree[k] for k in [ 
                "SimName", "Redshift", "H0", "Omega_M", "Omega_DE", "omega_b", 
                "w0", "wa", "n_s", "SODensity", "BoxSize", "ParticleMassHMsun", 
                "BoxSizeMpc", "ParticleMassMsun",
            ]
        } | {
            "sigma8" : AbacusCosmologySigma8Table[
                re.search(
                    r"AbacusSummit_\w+_c(\d\d\d)_ph\d\d\d", 
                    tree["SimName"], 
                ).group(1)
            ], 
            "NRange1" : list( mrange1 ), 
            "NRange2" : list( mrange2 ), 
        }, 
        "data" : {
            "r"  : rCenter, 
            "xi" : xi,
        }
    })
    af.write_to(output_file, all_array_compression = "zlib")
    af.close()

    logger.info("done!...")
    return

if __name__ == "__main__":
    corrcalc()
