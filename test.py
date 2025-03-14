import warnings
warnings.catch_warnings(action="ignore")

import asdf, glob
import numpy as np
import matplotlib.pyplot as plt
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

def testGalaxyCounts():

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

    centralCount   = np.zeros(massH.shape)
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
            _centralCount, _ = np.histogram( centralMass, bins = massBinEdges )
            centralCount[:] += _centralCount[:]
            
            # Number of satellite galaxies in each mass bin
            _satelliteCount, _ = np.histogram( centralMass[inverseIndex], bins = massBinEdges )
            satelliteCount[:] += _satelliteCount[:] - _centralCount[:]

    # -- 

    from abacus_galaxies import _loadCatalog

    files     = glob.glob( "../test-data/z3.000/halo_info/halo_info_*.asdf" )
    haloCount = np.zeros(massH.shape) 
    for file in files:
        # Halo mass in Msun
        _, _, haloMass = _loadCatalog(file)

        # Number of halos in each mass bin
        _haloCount, _ = np.histogram( haloMass * hubble, bins = massBinEdges )
        haloCount[:] += _haloCount[:]

    nonzeroMask    = (centralCount > 0.) & (haloCount > 0.)
    massH          = massH[nonzeroMask]
    satelliteCount = satelliteCount[nonzeroMask] / centralCount[nonzeroMask]
    centralCount   = centralCount[nonzeroMask] / haloCount[nonzeroMask]

    fig = plt.figure()
    plt.semilogx( massH, centralCount, 's', color = "tab:blue", label = "estimate Ncen" )
    plt.semilogx( massH, hm.centralCount( np.log(massH / hubble) ), '-.', color = "tab:blue", label = "expected Ncen" )
    plt.semilogx( massH, satelliteCount, 's', color = "red", label = "estimate Nsat" )
    plt.semilogx( massH, hm.satelliteCount( np.log(massH / hubble) ), '-.', color = "red", label = "expected Nsat" )
    plt.xlabel("Halo mass [Msun/h]")
    plt.ylabel("Satellite count")
    plt.xlim([1e+12, 1e+14])
    plt.ylim([0.0  , 4.0  ])
    plt.legend(ncols = 2)
    plt.grid()
    # plt.show()
    fig.savefig("../test-data/images/galaxy_count_hugebase_z3.png")

    return

def testSatellieMassFunction():

    files = glob.glob( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )

    with asdf.open(files[0]) as af:
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

    fracBinEdges = np.linspace(0.05, hm.scaleSHMF, 21)
    fracMass     = ( fracBinEdges[1:] + fracBinEdges[:-1] ) / 2.

    satelliteMassFunction = np.zeros(fracMass.shape)
    for file in files:
        with asdf.open(file) as af:

            uniqueIDs, centralIndex, inverseIndex = np.unique( 
                af["data"]["parentHaloID"], 
                return_index   = True, 
                return_inverse = True, 
            )

            # Mass of halos with central galaxy
            haloMass = af["data"]["galaxyMass"][centralIndex] 
            
            # Satellite galaxies
            satelliteIndex, = np.where( af["data"]["galaxyType"] == 2 )

            # Mass of the satllites as a fraction of the halo / central galaxy mass
            massFraction = af["data"]["galaxyMass"] / haloMass[inverseIndex]
            massFraction = massFraction[ massFraction < 1. ]
            
            _satelliteMassFunc, _     = np.histogram( massFraction, bins = fracBinEdges )
            satelliteMassFunction[:] += _satelliteMassFunc[:]

    satelliteMassFunction /= np.sum( satelliteMassFunction * np.diff(fracBinEdges) )

    # y = (fracMass / hm.scaleSHMF)**(-hm.slopeSHMF+1) * np.exp(-fracMass / hm.scaleSHMF) 

    plt.semilogy()
    plt.plot( fracMass, satelliteMassFunction, 's', color = "black" )
    # plt.plot( fracMass, y, '-.', color = "tab:blue" )
    plt.show()

    return

if __name__ == "__main__":
    # testGalaxyCatalogGeneration()
    # testGalaxyCounts()
    testSatellieMassFunction()



# a, b = 0.1, 0.5
# p = 1.91
# u = np.random.uniform(size = 10000)
# x = ( -( u*b**p - u*a**p - b**p ) / ( a*b )**p )**( -1./p )
# _, xb, _ = plt.hist(x, bins=21, density=True)
# xb = np.linspace(a, b, 51)
# yb = p * a**p * xb**(-p-1) / (1 - (a/b)**p)
# plt.plot(xb, yb, '--')
# plt.show()