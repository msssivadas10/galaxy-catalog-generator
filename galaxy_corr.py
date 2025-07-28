#!/usr/bin/python3
# 

__version__ = "0.1a"

import os, os.path, glob, logging, click, numpy as np, asdf
from typing import IO, Literal

@click.group
@click.version_option(__version__, message = "Galaxy correlation calculator - %(version)s")
def cli() -> None:
    r"""
    Tools for calculating galaxy correlation from abacus simulation based catalogs.
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
@click.option("--path", 
              type     = click.Path(file_okay = True, resolve_path = True, ), 
              required = True, 
              help     = "Path to the directory containing galaxy catalog files", )
@click.option("--mrange", 
              type     = ( click.FloatRange(0, min_open = True), )*2, 
              required = True, 
              help     = "Mass range for the parent halos (or central galaxies)", )

@click.option("--outfile", 
              type     = click.File("wb"), 
              required = True, 
              help     = "Path to the file to save new galaxy catalog", )
def savecat(
        path    : IO, 
        mrange  : tuple[float, float],
        outfile : IO, 
    ) -> None:
    r"""
    Load data from the catalog files and save a subset as a single file.
    """

    logger = logging.getLogger()

    if os.path.isdir(path):
        path = os.path.join( path, "galaxy_info_*.asdf" )
    files = glob.glob(path)
    if not files:
        return logger.warning("no catalog files in specified path")
    logging.info(f"found {len(files)} catalog files")
    
    header = {}
    with asdf.open(files[0]) as af:
        header.update( af["header"] )
        header["halosProcessed"]    = -1 # not needed - only for compatibility
        header["centralGalaxies"]   =  0
        header["satelliteGalaxies"] =  0
    header.update({ "NRange": list(mrange) })

    totalCount = 0
    ma, mb     = np.array(mrange) * header["particleMass"] # in Msun/h
    galaxyPositions, galaxyMass, parentHaloID, galaxyType = [], [], [], []
    for _idx, file in enumerate(files[:2], start = 1):
        logger.info(f"loading data from file {_idx} of {len(files)}: {file!r}")

        with asdf.open(file) as af:
            totalCount += af["header"]["centralGalaxies"]
            
            _parentHaloID    = af["data"]["parentHaloID"] 
            _galaxyPositions = af["data"]["galaxyPosition"] 
            _galaxyMass      = af["data"]["galaxyMass"] 
            _galaxyType      = af["data"]["galaxyType"] 

            start, = np.where( _galaxyType == 1 )
            stop   = np.hstack([ start[1:],  _galaxyType.shape[0] ])
            
            # Select halos in given mass range:
            haloMass    = _galaxyMass[start] # halo mass in Msun/h (same as central galaxy mass)
            filter_     = ( haloMass >= ma ) & ( haloMass <= mb ) 
            start, stop = start[filter_], stop[filter_]
            selection   = np.zeros(_galaxyType.shape, bool)
            for i, j in zip(start, stop):
                selection[i:j] = True

            _parentHaloID    = _parentHaloID[selection]
            _galaxyPositions = _galaxyPositions[selection]
            _galaxyMass      = _galaxyMass[selection]
            _galaxyType      = _galaxyType[selection]

            if _galaxyPositions.shape[0]:
                _, (nCentral, nSatellite) = np.unique(_galaxyType, return_counts = True)
                header["centralGalaxies"]   += nCentral
                header["satelliteGalaxies"] += nSatellite

                parentHaloID.append(_parentHaloID)
                galaxyPositions.append(_galaxyPositions)
                galaxyMass.append(_galaxyMass)
                galaxyType.append(_galaxyType)

    parentHaloID    = np.hstack(parentHaloID)
    galaxyPositions = np.vstack(galaxyPositions) # galaxy position in Mpc/h
    galaxyMass      = np.hstack(galaxyMass)      # galaxy mass in Msun/h
    galaxyType      = np.hstack(galaxyType)

    logger.info(f"saving catalog to file: {outfile.name!r}") 
    with asdf.AsdfFile({
            "header" : header, 
            "data" : {
                "parentHaloID"   : parentHaloID, 
                "galaxyPosition" : galaxyPositions, 
                "galaxyMass"     : galaxyMass, 
                "galaxyType"     : galaxyType,
            }
        }) as af:
        af.write_to(outfile, all_array_compression = 'zlib')
    
    return

@cli.command
@click.argument("cat1", type = click.File("rb"), )
@click.argument("cat2", type = click.File("rb"), required = False, )
@click.option("-o", "--outfile", 
              type     = click.File("wb"), 
              required = True, 
              help     = "Filename for the output (ASDF format)", )
@click.option("-r", "--rrange", 
              type     = ( click.FloatRange(0, min_open = True), )*2,  
              required = True,
              help     = "Distance bins range in Mpc/h", )
@click.option("--estimator", 
              type       = click.Choice(["nat", "ls", "ham", "dp"], case_sensitive = False), 
              default    = "nat", 
              help       = "Specify the estimator to use", )
@click.option("--rbins", 
              type     = click.IntRange(min = 2), 
              default  = 32, 
              help     = "Number of distance bins", )
@click.option("--nthreads", 
              type     = click.IntRange(min = 0), 
              default  = 0, 
              help     = "Number of threads to use (0 for using all available threads)", )
def galaxycorr(
        cat1      : IO, 
        cat2      : IO, 
        outfile   : IO, 
        rrange    : tuple[float, float],
        estimator : Literal["nat", "ls", "ham", "dp"] = "nat",
        rbins     : int = 16, 
        nthreads  : int = 0,
    ) -> None:
    r"""
    Calculate galaxy correlation between given catalogs.

    NOTE: This assumes both catalogs are derived from a single catalog - that is, model 
    parameters and boxsize, redshift etc should be same.

    """
    logger = logging.getLogger()

    logger.info(f"loading catalog 1, {cat1.name!r}")
    with asdf.open(cat1) as af:
        tree    = af["header"] 
        mrange1 = af["header"]["NRange"]
        D1 = af["data"]["galaxyPosition"].copy()
        T1 = af["data"]["galaxyType"].copy()

    D2, T2, mrange2 = None, None, mrange1
    if cat2 is not None:
        logger.info(f"loading catalog 2, {cat2.name!r}")
        with asdf.open(cat2) as af: 
            mrange2 = af["header"]["NRange"]
            D2 = af["data"]["galaxyPosition"].copy()
            T2 = af["data"]["galaxyType"].copy()

    # -- Calculating correlation function

    from galaxy_catalog_generator.misc.correlation import PairCountData

    estimator = estimator.lower()
    rBinEdges = np.logspace( np.log10(rrange[0]), np.log10(rrange[1]), rbins+1 ) # in Mpc/h
    rCenter   = np.sqrt( rBinEdges[1:] * rBinEdges[:-1] ) 

    # Galaxy-galaxy correlation
    countData = PairCountData.countPairs(
        D1, 
        D2, 
        rBinEdges, # in Mpc/h
        boxsize  = tree["boxsize"], 
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
    logger.info(f"saving correlation values to file {outfile.name!r}")
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
        af.write_to(outfile, all_array_compression = "zlib")

    logger.info("correlation calculation completed!...")
    return

if __name__ == "__main__":
    cli()