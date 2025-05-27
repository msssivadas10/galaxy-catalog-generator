#!/usr/bin/python3
# 
# CLI tools to process abacus halo catalogs and generate object count grids.
# 

__version__ = "0.1a"

ARRAY_SIZE_LIMIT = 536870912 # maximum size limit for count file: 512 MiB

import os, os.path, glob, re, click, logging, numpy as np
from typing import TextIO, Literal, Any
from numpy.typing import NDArray
from abacusnbody.data.read_abacus import read_asdf
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog

############################################################################################################
#                                   HELPER FUNCTIONS FOR CALCULATIONS                                      #
############################################################################################################

def listCatalogs(simname: str, redshift: float, particle: Literal["halos", "field"], search_paths: list[str] = []) -> tuple[list[str], Literal["halos", "field"]]:
    r"""
    Search for abacus particle catalog files for given simulation name and redshift in the paths
    specified and return a list of paths to the files.  
    """

    logger = logging.getLogger()
    
    # Going through all search paths, looking for catalog files, until files are found...
    # If no search path specified, look in the current directory.
    prefix = dict( halos = "halo_info", field = "field_rv_A" ).get(particle, None)
    if not prefix:
        logger.error(f"incorrect particle type {particle!r}")
        return [], "unknown"
    tailPath   = os.path.join(simname, "halos", f"z{redshift:.3f}", prefix , f"{prefix}_*.asdf")
    globResult = []
    for path in search_paths:
        globResult = glob.glob( os.path.join(path, tailPath) )
        if globResult: break

    # Filter out only the files with names matchng the pattern and sorting based on the
    # integer index: 
    _files = { 
        int( m.group(1) ): m.string 
            for m in map( 
                lambda fn: re.search(prefix + r"_(\d{3}).asdf", fn), # for files like `halo_info_123.asdf`
                globResult 
            ) 
            if m 
    }
    files  = [ _files[idx] for idx in sorted(_files) ]

    if not files: 
        logger.info(f"no halo files for simulation {simname!r} at redshift {redshift}: exiting...")
    else:
        logger.info(f"found {len(files)} halo files for simulation {simname!r} at redshift {redshift}")

    return files, particle

def loadCatalog(fn: str, particle: Literal["halos", "field"]) -> NDArray[np.float64]:
    r"""
    Load data from a abacus particle catalog file. Particle position in Mpc/h and mass in particle unit are 
    returned. Field particles have no mass output.
    """

    logger = logging.getLogger()
    
    logger.info(f"loading halo catalog from file: {fn!r}")
    if particle == "halos":
        # Load halo catalog
        catalog = CompaSOHaloCatalog(fn, cleaned = False, fields = ["id", "SO_central_particle", "N"])
        data    = np.hstack([
            # -- Halo position coordinates in Mpc/h units
            np.array( catalog.halos["SO_central_particle"] ),
            # -- Halo mass as particle count
            np.array( catalog.halos["N"], dtype = np.float64 ).reshape(-1, 1)
        ])
        
    else:
        # Load field catalog
        # -- Particle position coordinates in Mpc/h units
        data = np.array( read_asdf(fn, [ "pos" ])[ "pos" ] )


    # Sorting the objects based on x coordinate
    sortedOrder = np.argsort( data[:,0] )
    data        = data[sortedOrder, :]

    return data

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

############################################################################################################
#                                               CLI AND COMMANDS                                           #
############################################################################################################

@click.command
@click.version_option(__version__, message = "Abacus Count %(version)s")
@click.option("--siminfo", 
              type     = (str, float), 
              required = True, 
              callback = siminfoValidator, 
              help     = "A tuple of a valid abacus simulation name and redshift value", )
@click.option("-s", "--cell-size", 
              type     = click.FloatRange(min = 0., max = 1., min_open = True), 
              required = True, 
              help     = "Size of the cell as a fraction of the boxsize", )
@click.option("--output-path", 
              type     = click.Path(file_okay = False), 
              required = True, 
              help     = "Path to the directory to save count data files", )
@click.option("-g", "--cell-geometry", 
              type       = click.Choice(["square", "sphere"], case_sensitive = False), 
              required = True, 
              help       = "Specify the shape of the cells to use", )
@click.option("-t", "--particle", 
              type       = click.Choice(["halos", "field"], case_sensitive = False), 
              required = True, 
              help       = "Specify the particle to use", )
@click.option("-b", "--massbins", 
              type     = click.File("r"), 
              help     = "Name of the file containing mass bin edges", )
@click.option("-p", "--path" , 
              type     = click.Path(exists = True), 
              default  = [ os.getcwd() ], 
              multiple = True,  
              help     = "Paths to look for simulation catalog files", )
def count(
        siminfo       : tuple[str, float], 
        cell_size     : float, 
        output_path   : str,
        cell_geometry : Literal["square", "sphere"] = "square",
        particle      : Literal["halos", "field"]   = "halos",
        massbins      : TextIO                      =  None,
        path          : list[str]                   =  [], 
    ) -> None:
    r"""
    Calculate count-in-cells.

    Calculate count-in-cells from abacus summit simulation catalogs. Count both halo and field particles.
    """

    import logging.config, warnings, asdf
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
    
    logger = logging.getLogger()

    # List all catalog files in the given path:
    simname, redshift = siminfo
    files, particle = listCatalogs(simname, redshift, particle, path)
    if not files: return
    
    
    # If the output path does not exist, create...
    output_path = os.path.join(output_path, simname, f"z{redshift:.3f}", f"{particle}_count")
    os.makedirs(output_path, exist_ok = True)

    # Load mass binedges from the given file
    massBinEdges, massBins = np.array([ 0., np.inf ]), 1
    if massbins and particle == "halos":
        with massbins:
            massBinEdges = np.loadtxt(massbins, comments = '#')
            massBins     = len(massBinEdges) - 1

    # Header: contains the values of cosmology and simulation parameters
    import abacusnbody.metadata

    logger.info(f"loading metatdata for simulation {simname!r}for redshift {redshift}...")
    header = abacusnbody.metadata.get_meta(simname, redshift)
    header = { 
        _key: header[_key] for _key in [
            "Redshift", "SimName", "H0", "Omega_M", "Omega_DE", "omega_b", "w0", "wa", "n_s",
            "SODensity", "BoxSize", "ParticleMassHMsun", "BoxSizeMpc", "ParticleMassMsun", 
    ]}
    header["CellsizeFraction"] = cell_size
    
    # Calculate the actual cellsize and shape of count array
    boxsize  = header["BoxSize"]
    cellsize = cell_size * boxsize
    header["CellSize"]    = cellsize
    header["CellSizeMpc"] = cell_size * header["BoxSizeMpc"]
    header["CellShape"]   = cell_geometry

    gridsize = int( boxsize / cellsize )
    logger.info(f"using cellsize {cellsize:3g} Mpc/h and {gridsize}^3 cells...")

    # Based on the total size of the array, the cell array is partitioned into multiple slabs along
    # X direction, so that the total size of a slab will always be < 500 MiB. This will make sure the 
    # enough memory available, as most systems will have at least 4 GB of RAM :)   
    totalArraySize   = gridsize**3 * massBins * 8 # total array size in bytes
    if totalArraySize > ARRAY_SIZE_LIMIT:
        # Partition the array into multiple slabs along x axis - find the new x axis grid size to 
        # satisfy the condition...
        slabCount     = totalArraySize // ARRAY_SIZE_LIMIT # no. of slabs required
        newGridsizeX  = gridsize // slabCount
        remGridsizeX  = gridsize - newGridsizeX * slabCount
    else:
        # No partition is required
        slabCount     = 1
        newGridsizeX  = gridsize
        remGridsizeX  = 0 
    
    # Cell edges: cells along y and z axis are the same, and span the entire space. Cells
    # along the x axis is different for different partitions...  
    qedges     = np.arange(gridsize+1) * cellsize - 0.5*boxsize # cell edges along y or z
    xedgesList = [] # cell edges along x axis

    # Create the partitions and initialize the count to 0, save as ASDF:
    partitionSize, partitionEdges, outputFiles = [], [ -0.5*boxsize ], []
    for i in range(slabCount):
        _gridsize = newGridsizeX
        if remGridsizeX > 0:
            _gridsize    += 1
            remGridsizeX -= 1

        # X range for the partition:
        xa = partitionEdges[-1]
        xb = xa + _gridsize * cellsize
        
        partitionSize.append(_gridsize)
        partitionEdges.append(xb)

        xedges = np.arange(_gridsize+1) * cellsize + xa 
        xedgesList.append(xedges)

        with asdf.AsdfFile({
                "header" : {
                    **header, 
                    "PartitionIndex" : i, 
                    "PartitionRange" : [ xa, xb ],
                }, 
                "data" : {
                    "count"  : np.zeros(( _gridsize, gridsize, gridsize, massBins ), dtype = np.int64), 
                    "xedges" : xedges, 
                    "yedges" : qedges,
                    "medges" : massBinEdges, # NOTE: mass bins only makes sense for halo particles...
                },
            }) as af:
            file = os.path.join(output_path, f"{particle}_count_{i:03d}.asdf")
            af.write_to( file )
            outputFiles.append( file )

    # ====================================== Count Calculation =========================================== #

    for fn in files:

        # Load position and mass (for halos only)...
        data   = loadCatalog(fn, particle) # [ ...position, mass ]
        xa, xb = data[0, 0], data[-1, 0]   # range of x values in the position array
        
        for i, outputFile in enumerate(outputFiles): 

            if not ( xa <= partitionEdges[i+1] and partitionEdges[i] <= xb ):
                continue # range is outside actual x coordinate range

            # Filter for spherical cell: remove points outside the sphere
            # The spherical cells are actually inside the square cells, which makes the counting easy, 
            # as we only need to remove points that are far away than cellsize/2 distance from the 
            # cell center... 
            if cell_geometry == "sphere":
                logger.info(f"filtering out objects outside spherical cells...")

                _data = ( data[:, :3] + 0.5*boxsize ) / cellsize
                _data = _data - np.floor(_data) - 0.5  # position relative to the cell center
                _data = np.sum( _data**2, axis = 1 )   # squared distance from the cell center
                
                # Filtered data, containing only points inside circular cells
                _data = data[ _data < 0.25 ]
            else: # square
                _data = data

            logger.info(f"counting in partition #{i}...")
            if particle == "field": # field particle (no mass column)
                arr, _ = np.histogramdd(
                    _data, 
                    bins    = [ xedgesList[i], qedges, qedges ], 
                    density = False, 
                ) 
                arr = np.expand_dims(arr, axis = -1)
            else: # halo particles
                arr, _ = np.histogramdd(
                    _data, 
                    bins    = [ xedgesList[i], qedges, qedges, massBinEdges ], 
                    density = False, 
                )
            arr = arr.astype( np.int64 )

            logger.info(f"updating data in {outputFile!r}...")
            with asdf.open( outputFile, mode="rw" ) as af:
                filearr     = af['data']["count"]
                filearr[:] += arr[:] 
                af.update() # save changes 
        
    logger.info("count calculations completed! :)")
    return

if __name__ == "__main__":
    count()
