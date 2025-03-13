import numpy as np
from astropy.cosmology import FlatLambdaCDM, w0waCDM
from galaxy_catalog_generator.halomodel import HaloModel

def testGalaxyCatalogGeneration():
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

def testGalaxyCatalog():
    import asdf, glob
    import matplotlib.pyplot as plt

    files = glob.glob( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )

    with asdf.open(files[0]) as af:
        unitMass = af["header"]["ParticleMassHMsun"] 
        hubble   = 0.01 * af["header"]["H0"]

        # Halo model object
        hm = HaloModel.create(
            Mmin      = af["header"]["HaloModel_Mmin"     ], 
            sigmaM    = af["header"]["HaloModel_sigmaM"   ], 
            M0        = af["header"]["HaloModel_M0"       ], 
            M1        = af["header"]["HaloModel_M1"       ], 
            alpha     = af["header"]["HaloModel_alpha"    ], 
            scaleSHMF = af["header"]["HaloModel_scaleSHMF"], 
            slopeSHMF = af["header"]["HaloModel_slopeSHMF"], 
            redshift  = af["header"]["Redshift"],
            cosmo     = w0waCDM(
                H0    = af["header"]["H0"], 
                Om0   = af["header"]["Omega_M"], 
                Ode0  = af["header"]["Omega_DE"],
                Ob0   = af["header"]["omega_b"] / hubble**2, 
                w0    = af["header"]["w0"], 
                wa    = af["header"]["wa"],
                Tcmb0 = af["header"]["CMBTemperature"], 
            ), 
            psmodel   = "eisenstein98_zb", 
            mfmodel   = "tinker08", 
            ns        = af["header"]["n_s"], 
            sigma8    = af["header"]["sigma8_m"], 
            Delta     = af["header"]["SODensity"][0],
        )

    massBinEdges = unitMass * np.logspace(1.2, 3, 21) 
    massH        = np.sqrt( massBinEdges[1:] * massBinEdges[:-1] )

    haloCount      = np.zeros(massH.shape)
    satelliteCount = np.zeros(massH.shape)
    for file in files:
        with asdf.open(file) as af:

            uniqueIDs, centralIndex, inverseIndex = np.unique( 
                af["data"]["parentHaloID"], 
                return_index   = True, 
                return_inverse = True, 
            )

            # Mass of halos with central galaxy
            centralMass = af["data"]["galaxyMass"][centralIndex] 
            
            # Number of halos with central galaxies in each mass bin
            _haloCount, _ = np.histogram( centralMass, bins = massBinEdges )
            haloCount[:] += _haloCount[:]
            
            # Number of satellite galaxies in each mass bin
            _satelliteCount, _ = np.histogram( centralMass[inverseIndex], bins = massBinEdges )
            satelliteCount[:] += _satelliteCount[:] - _haloCount[:]


    nonzeroMask    = (haloCount > 0.)
    massH          = massH[nonzeroMask]
    satelliteCount = satelliteCount[nonzeroMask] / haloCount[nonzeroMask]

    plt.semilogx( massH, satelliteCount, 's', color = "black", label = "estimate" )
    plt.semilogx( massH, hm.satelliteCount( np.log(massH / hubble) ), '-', color = "tab:blue", label = "expected" )
    plt.xlabel("Halo mass [Msun/h]")
    plt.ylabel("Satellite count")
    plt.xlim([1e+12, 1e+14])
    plt.legend()
    plt.show()

    return

if __name__ == "__main__":
    testGalaxyCatalog()
    # testGalaxyCatalogGeneration()
