
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from galaxy_catalog_generator.halomodel import HaloModel

def main():

    hm = HaloModel.create(
        Mmin      = 1e+08, 
        sigmaM    = 0., 
        M0        = 1e+06, 
        M1        = 1e+12, 
        alpha     = 0.3, 
        scaleSHMF = 0.39, 
        slopeSHMF = 1.91, 
        redshift  = 0.,
        cosmo     = FlatLambdaCDM(H0 = 70., Om0 = 0.3, Ob0 = 0.05, Tcmb0 = 2.725), 
        psmodel   = "eisenstein98_zb", 
        mfmodel   = "tinker08", 
        ns        = 1., 
        sigma8    = 1., 
        Delta     = 200,
    )

    haloMass = 1e+13 # in Msun
    satelliteCatalog = hm.generateSatellitePositions( np.log(haloMass), [10., 10., 10.] )
    print(satelliteCatalog)
    return

if __name__ == "__main__":
    main()

# NOTE: Generating a complete satellite galaxy catalog from a halo catalog
# 
# haloCatalog = # load the catalog file...
# galaxyCatalog = [] 
# for haloPosition, haloMass in haloCatalog:
#   satelliteCatalog = hm.generateSatellitePositions( np.log(haloMass) )
#   satelliteCatalog[:,0:3] += haloPosition 
#   galaxyCatalog.append(satelliteCatalog)
#   # galaxyCatalog.append([[ central galaxy X, Y, Z and mass ]]) 
# galaxyCatalog = np.concatenate(galaxyCatalog, axis = 0)
# # save galaxyCatalog