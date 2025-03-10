
import os, os.path
# import logging, logging.config
# import asdf, click, yaml
import numpy as np
from numpy.typing import NDArray
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from galaxy_catalog_generator.halomodel import HaloModel

def _listCatalogs(path: str, /) -> list[str]:
    r"""
    Return a list of halo catalog files for the given path or glob pattern.
    """
    
    import glob, re

    globResult = glob.glob(path)
    if not globResult:
        return []
    
    filePattern = re.compile(r"halo_info_(\d+).asdf") # for files like `halo_info_123.asdf`
    files       = [
        file for file, _ in sorted(
            # Filter out only the files with names matchng the pattern and sorting based on the
            # integer index: 
            filter(
                lambda a: a[1], 
                [ 
                    (file, filePattern.match(os.path.basename(file))) for file in globResult 
                ]
            ), 
            key  = lambda a: int( a[1].group(1) ), 
            )
    ]
    return files

def _loadCatalog(fn: str, /) -> NDArray[np.float64]:
    r"""
    Read a halo catalog file and return the halo position (in Mpc/h) and mass (in Msun/h). Position 
    coordinates are the first 3 columns of the array and mass is the last column. 
    """
    
    catalog = CompaSOHaloCatalog(fn, cleaned = False, fields = ["SO_central_particle", "N"])
    
    # Halo position coordinates in Mpch units
    haloPosition = np.array( catalog.halos["SO_central_particle"] )
    
    # Halo mass in Msun/h
    unitMass = catalog.header["ParticleMassHMsun"]
    haloMass = np.array( catalog.halos["N"] ) * unitMass
    
    data = np.concatenate([ haloPosition, haloMass[:,None] ], axis = 1)
    return data

def generateGalaxyCatalog(
        path     : str, 
        # Mmin     : float,
        # sigmaM   : float,
        # M0       : float,
        # M1       : float,
        # alpha    : float,
        # scaleSHMF: float,
        # slopeSHMF: float,
    ) -> HaloModel:

    files = _listCatalogs(path)
    if not files: return

    import asdf

    # Create a halo model using the cosmology used for simulation (loaded from the header section of
    # the first file) and other given parameters...
    header = dict.fromkeys([
        "SODensity", "BoxSize", "H0", "Omega_DE", "Omega_K", "Omega_M", "omega_b", "n_s", "w0", "wa", 
        "Redshift", "SimName", 
    ])
    with asdf.open(files[0]) as af:
        for key, value in dict.items( af["header"] ):
            if key not in header:
                continue
            header[key] = value

    # haloModel = HaloModel.create(
    #     Mmin, 
    #     sigmaM, 
    #     M0, 
    #     M1, 
    #     alpha, 
    #     scaleSHMF, 
    #     slopeSHMF, 
    #     redshift, 
    #     cosmo, 
    #     psmodel = "eisenstein98_zb", 
    #     mfmodel = "tinker08", 
    #     ns = header["n_s"], 
    #     sigma8 = ..., 
    #     Delta = ...,
    # )
    return


path  = "/home/ms3/Documents/phd/cosmo/workspace/cosmo-calc/workspace/z2.000/halo_info/halo_info_*.asdf"
generateGalaxyCatalog(path)