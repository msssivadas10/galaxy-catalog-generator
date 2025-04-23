
import numpy as np
import click, yaml, inspect, logging.config

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
        hubble   = 0.01 * af["header"]["H0"]
        boxsize  = af["header"]["BoxSizeMpc"]

    # Estimating halo mass-function
    massBinEdges = unitMass * np.logspace( np.log10(start), np.log10(stop), bins ) 
    dlnMassH     = np.diff( np.log( massBinEdges ) )
    dndM         = np.zeros(dlnMassH.shape)
    for file in files:
        logger.info(f"loading halo catalog from file: {file!r}")
        _, _, haloMass = _loadCatalog(file) # halo mass in Msun 
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