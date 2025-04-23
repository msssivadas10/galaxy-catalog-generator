__version__ = "0.1a"


import logging, logging.config
import os, os.path, glob, re, asdf, click
import numpy as np
from numpy.typing import NDArray
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog

def _listCatalogs(path: str, /) -> list[str]:
    r"""
    List all the files matching to the given path or glob pattern.
    """
    
    globResult  = glob.glob(path)
    filePattern = re.compile(r"halo_info_(\d+).asdf") # for files like `halo_info_123.asdf`
    files       = [
        file for file, _ in sorted(
            # Filter out only the files with names matchng the pattern and sorting based on the
            # integer index: 
            filter(
                lambda a: a[1], 
                [ ( file, filePattern.match(os.path.basename(file)) ) for file in globResult ]
            ), 
            key  = lambda a: int( a[1].group(1) ), 
        )
    ]
    return files

def _loadHaloMass(fn: str, /) -> NDArray[np.float64]:
    r"""
    Load an abacus summit halo catalog file and return the halo mass (Msun).
    """
    
    catalog  = CompaSOHaloCatalog(fn, cleaned = False, fields = ["N"])
    
    # Halo mass in Msun
    unitMass = catalog.header["ParticleMassMsun"]
    haloMass = np.array( catalog.halos["N"] ) * unitMass

    return haloMass

@click.command
@click.argument('path' , type = click.Path(exists = False))
@click.option("--out"  , type = click.Path(exists = False), default = "hmf.csv")
@click.option("--start", type = float, default = 2     , help = "Minimum mass as number of particle")
@click.option("--stop" , type = float, default = 10_000, help = "Maximum mass as number of particle")
@click.option("--bins" , type = int  , default = 32    , help = "Number of mass bins"               )
@click.option("--warn" , is_flag = False, help = "Enable or disable warnings")
def hmf(
        path : str, 
        out  : str,
        start: float,
        stop : float,
        bins : int,
        warn : bool = False,
    ) -> None:
    r"""
    Calculate halo mass-function from abacus summit halo catalogs.
    """

    if not warn:
        import warnings
        warnings.catch_warnings(action="ignore")

    # Configure logger
    logging.config.dictConfig({
            'version': 1, 
            'disable_existing_loggers': True, 
            'formatters': {
                'default': { 'format': '[ %(asctime)s %(levelname)s %(process)d ] %(message)s' }
            }, 
            'handlers': {
                'stream': {
                    'level': 'INFO', 
                    'formatter': 'default', 
                    'class': 'logging.StreamHandler', 
                    'stream': 'ext://sys.stdout'
                }, 
                'file': {
                    'level': 'INFO', 
                    'formatter': 'default', 
                    'class': 'logging.handlers.RotatingFileHandler', 
                    'filename': '.'.join([ os.path.basename(__file__).rsplit('.', 1)[0], "log" ]), 
                    'mode': 'a', 
                    'maxBytes': 2097152, # create a new file if size exceeds 2 MiB
                    'backupCount': 4     # use maximum 4 files
                }
            }, 
            'loggers': {
                'root': {
                    'level': 'INFO', 
                    'handlers': [
                        'stream', 
                        'file'
                    ]
                }
            }
        })
    logger = logging.getLogger()

    # List all catalog files in the given path:
    files = _listCatalogs(path)
    logger.info(f"listed {len(files)} files with path pattern {path!r}")
    if not files: return

    # Read the first catalog file and load the header data. 
    logger.info(f"loading header data from file: {files[0]!r}...")
    with asdf.open(files[0]) as af:
        unitMass = af["header"]["ParticleMassMsun"] 
        boxsize  = af["header"]["BoxSizeMpc"]

    # Estimating halo mass-function
    massBinEdges = unitMass * np.logspace( np.log10(start), np.log10(stop), bins ) 
    dlnMassH     = np.diff( np.log( massBinEdges ) )
    dndM         = np.zeros(dlnMassH.shape)
    for file in files:
        logger.info(f"loading halo catalog from file: {file!r}")
        haloMass  = _loadHaloMass(file) # halo mass in Msun 
        dndM[:]  += np.histogram(haloMass, bins = massBinEdges)[0][:]
    dndM  /= ( boxsize**3 * dlnMassH )
    massH = np.sqrt( massBinEdges[1:] * massBinEdges[:-1] )

    # Save data
    logger.info(f"saving data to file: {out!r}")
    np.savetxt(
        out, 
        np.log( np.stack((massH, dndM), axis = -1) ), 
        header    = "mass (Msun), dndm (Mpc^-3)",
        delimiter = ',',
        comments  = '#'
    )
    return

if __name__ == "__main__":
    hmf()