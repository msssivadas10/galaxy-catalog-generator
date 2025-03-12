
import os, os.path, glob, re, asdf
import numpy as np
from numpy.typing import NDArray
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
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

def listCatalogs(path: str, /) -> list[str]:
    
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

def loadCatalog(fn: str, /) -> tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64]]:
    
    catalog  = CompaSOHaloCatalog(fn, cleaned = False, fields = ["id", "SO_central_particle", "N"])
    
    # Halo position coordinates in Mpc units
    hubble       = 0.01 * catalog.header["H0"]
    haloPosition = np.array( catalog.halos["SO_central_particle"] ) / hubble
    
    # Halo mass in Msun
    unitMass = catalog.header["ParticleMassMsun"]
    haloMass = np.array( catalog.halos["N"] ) * unitMass

    # Halo integer id
    haloID = np.array( catalog.halos["id"] )
    
    return haloID, haloPosition, haloMass

def generateGalaxyCatalog(
        catpath      : str,
        outpath      : str,  
        Mmin         : float, 
        sigmaM       : float, 
        M0           : float, 
        M1           : float, 
        alpha        : float, 
        scaleSHMF    : float, 
        slopeSHMF    : float, 
        powerspectrum: str  = "eisenstein98_zb", 
        massfunction : str  = "tinker08", 
    ) -> None:

    # Check if the output path is existing or not: create it if not exist: NOTE: use paths to folder 
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # List all catalog files in the given path:
    files = listCatalogs(catpath)
    if not files: return

    # Read the first catalog file and load the header data. This is used for initialising cosmology 
    # nad halo model objects and is saved on each galaxy catalog files:
    catalogHeader = {}
    with asdf.open(files[0]) as af:
        catalogHeader.update( af["header"] )

        # Sigma-8 value is taken from the look-up table, using the cosmology ID
        catalogHeader["sigma8_m"] = _AbacusCosmologySigma8Table.get(
            re.search(
                r"AbacusSummit_\w+_c(\d\d\d)_ph\d\d\d",  
                catalogHeader["SimName"], 
            ).group(1), 
            None, 
        )

        # Other halo model parameters and models:
        catalogHeader["HaloModel_Mmin"     ] = Mmin
        catalogHeader["HaloModel_sigmaM"   ] = sigmaM
        catalogHeader["HaloModel_M0"       ] = M0
        catalogHeader["HaloModel_M1"       ] = M1
        catalogHeader["HaloModel_alpha"    ] = alpha
        catalogHeader["HaloModel_scaleSHMF"] = scaleSHMF
        catalogHeader["HaloModel_slopeSHMF"] = slopeSHMF
        catalogHeader["PowerSpectrumModel" ] = powerspectrum
        catalogHeader["MassFunctionModel"  ] = massfunction
        catalogHeader["CMBTemperature"     ] = 2.728
    
    # Using the parameters from the simulation metadata and given halo model parameters, create
    # the halo model object for generating galaxy catalogs:  
    haloModel = HaloModel.create(
        Mmin, 
        sigmaM, 
        M0, 
        M1,
        alpha, 
        scaleSHMF, 
        slopeSHMF, 
        redshift  = catalogHeader["Redshift"], 
        cosmo     = w0waCDM(
            H0    = catalogHeader["H0"], 
            Om0   = catalogHeader["Omega_M"], 
            Ode0  = catalogHeader["Omega_DE"],
            Ob0   = catalogHeader["omega_b"] / ( 0.01*catalogHeader["H0"] )**2,
            w0    = catalogHeader["w0"], 
            wa    = catalogHeader["wa"], 
            Tcmb0 = catalogHeader["CMBTemperature"],
        ), 
        psmodel = powerspectrum, 
        mfmodel = massfunction, 
        ns      = catalogHeader["n_s"], 
        sigma8  = catalogHeader["sigma8_m"], 
        Delta   = catalogHeader["SODensity"][0],
    )

    filePattern = re.compile(r"halo_info_(\d+).asdf") # for files like `halo_info_123.asdf`
    for file in files:
        parentHaloID, galaxyPositions, galaxyMass = [], [], []

        haloID, haloPosition, haloMass = loadCatalog(file)
        for id, posH, massH in enumerate(haloID, haloPosition, haloMass):
            galaxyData = haloModel.generateSatellitePositions( np.log(massH), posH )
            if galaxyData is None:
                continue

            parentHaloID.append( np.repeat(id, repeats = galaxyData.shape[0]) )
            galaxyPositions.append( galaxyData[:, 0:3] )
            galaxyMass.append( galaxyData[:, 3] )

        parentHaloID    = np.hstack(parentHaloID)
        galaxyPositions = np.vstack(galaxyPositions)
        galaxyMass      = np.hstack(galaxyMass)

        outfile = os.path.join(
            outpath, 
            f"galaxy_info_{ filePattern.match(os.path.basename(file)).group(1) }.asdf", 
        )
        af = asdf.AsdfFile({
            "header": catalogHeader, 
            "data"  : {
                "parentHaloID": parentHaloID, 
                "positionMpc" : galaxyPositions, 
                "massMsun"    : galaxyMass, 
            }
        })
        af.write_to(outfile, all_array_compression = 'zlib')
        af.close()
        
    return


# TODO: main function: cli


# XXX: **test**
generateGalaxyCatalog(
    catpath   = os.path.abspath("../cosmo-calc/workspace/z2.000/halo_info/halo_info_*.asdf"), 
    outpath   = os.path.abspath("../galaxy_catalog/z2.000/"), 
    Mmin      = 1e+08, 
    sigmaM    = 0., 
    M0        = 1e+06, 
    M1        = 1e+12, 
    alpha     = 0.3, 
    scaleSHMF = 0.39, 
    slopeSHMF = 1.91, 
)
