import warnings
warnings.catch_warnings(action="ignore")

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from galaxy_catalog_generator.halomodel import HaloModel
from test_catalog import setupLoading

def testGalaxyCatalogGeneration():
    # Test galaxy catalog generation

    hm = HaloModel.create(
        Mmin      = 1e+12, 
        sigmaM    = 0., 
        M0        = 1e+12, 
        M1        = 1e+12, 
        alpha     = 0.3, 
        scaleSHMF = 0.50, 
        slopeSHMF = 2.00, 
        redshift  = 3.,
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

def testMinimizeHaloModel():
    observedGalaxyDensity = 6.7e-04
    # observedGalaxyDensity = 2.5e-04
    observedSatelliteFraction = 0.1

    # Load halo model paameters from the galaxy catalog
    files, unitMass, hubble, hm = setupLoading( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )
    object.__setattr__(hm, "alpha", 1.)

    def scoreFunction(m):
        Mmin, M1 = np.exp(m[0]), np.exp(m[1])
        object.__setattr__(hm, "Mmin", Mmin)
        object.__setattr__(hm, "M0"  , Mmin)
        object.__setattr__(hm, "M1"  , M1  )
        galaxyDensity     = hm.averageGalaxyDensity( lnmb = np.log(1e+20) )
        satelliteFraction = hm.averageSatelliteFraction( lnmb = np.log(1e+20) )
        score1 = (galaxyDensity - observedGalaxyDensity)**2 / observedGalaxyDensity**2 
        score2 = (satelliteFraction - observedSatelliteFraction)**2 / observedSatelliteFraction**2
        return np.log(score1 + score2 + 1e-16)

    print("Making grid...")
    Mmin = np.logspace(12, 14, 21)
    M1   = np.logspace(12, 16, 21)
    Mmin, M1 = np.meshgrid(Mmin, M1)
    score = np.zeros_like(Mmin)
    for i in range( score.shape[0] ):
        for j in range( score.shape[1] ):
            score[i, j] = scoreFunction( np.log([ Mmin[i, j], M1[i, j] ]) )

    # Minimize to find best fitting HOD parameters
    i, j = np.unravel_index(np.argmin(score), score.shape)
    optMmin = Mmin[i, j]
    optM1   = M1[i, j]
    print(f"Minimum point (guess): Mmin = {optMmin:.2e} Msun, M1 = {optM1:.2e} Msun")
    
    object.__setattr__(hm, "Mmin", optMmin)
    object.__setattr__(hm, "M0"  , optMmin)
    object.__setattr__(hm, "M1"  , optM1  )
    print(f" - Galaxy density = {hm.averageGalaxyDensity():.3e}")
    print(f" - Satellite fraction = {hm.averageSatelliteFraction():.3g}")

    # print("Starting optimization for (Mmin, M1)...")
    from scipy.optimize import minimize

    optMmin = 5e+12 
    optM1   = 1e+14
    res = minimize(
        scoreFunction, 
        x0 = [ np.log(optMmin), np.log(optM1) ], 
        bounds = [( np.log(1e+12), np.log(1e+13) ), ( np.log(1e+13), np.log(1e+15) )]
    )
    # print( res )
    optMmin, optM1 = np.exp(res.x) 
    print(f"Minimum point (optimum): Mmin = {optMmin:.2e} Msun, M1 = {optM1:.2e} Msun")
    print(f"Status: success={res.success}")

    object.__setattr__(hm, "Mmin", optMmin)
    object.__setattr__(hm, "M0"  , optMmin)
    object.__setattr__(hm, "M1"  , optM1  )
    print(f" - Galaxy density = {hm.averageGalaxyDensity():.3e}")
    print(f" - Satellite fraction = {hm.averageSatelliteFraction():.3g}")

    fig = plt.figure()
    plt.loglog()
    plt.contourf(Mmin, M1, score, levels = 21)
    plt.plot([optMmin], [optM1], 'o', color = "black")
    plt.text(optMmin, optM1, f"  ({optMmin:.2e}, {optM1:.2e})")
    plt.colorbar()
    plt.xlabel("Mmin [Msun]")
    plt.ylabel("M1 [Msun]")
    fig.savefig("../test-data/images/mmin_m1_score.png")
    # plt.show()
    return

def testSatelliteFraction():

    files, unitMass, hubble, hm = setupLoading( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )
    object.__setattr__(hm, "alpha", 1.)

    MminList = np.logspace(8, 16, 21)
    M1List   = [ 1e+12, 1e+14, 1e+16 ]
    satelliteFractions = []
    for M1 in M1List:
        satelliteFraction = []
        for Mmin in MminList:
            object.__setattr__(hm, "Mmin", Mmin)
            object.__setattr__(hm, "M0"  , Mmin)
            object.__setattr__(hm, "M1"  , M1  )
            satelliteFraction.append( hm.averageSatelliteFraction() )
        satelliteFractions.append( np.array(satelliteFraction) )

    fig = plt.figure()
    plt.loglog()
    for satelliteFraction, M1 in zip(satelliteFractions, M1List):
        plt.plot(MminList, satelliteFraction, "-.", zorder = 2, label = f"{M1:.2e}")
    plt.axhline(0.1, color = "black", zorder = 1)
    plt.legend(title = "M1 [Msun]", ncols = 3)
    plt.xlabel("Mmin [Msun]")
    plt.ylabel("fsat")
    plt.grid()
    fig.savefig("../test-data/images/satellite_fraction.png")
    # plt.show()
    return

def testGalaxydensity():

    files, unitMass, hubble, hm = setupLoading( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )
    object.__setattr__(hm, "alpha", 1.)

    MminList = np.logspace(8, 16, 21)
    M1List   = [ 1e+12, 1e+14, 1e+16 ]
    galaxyDensities = []
    for M1 in M1List:
        galaxyDensity = []
        for Mmin in MminList:
            object.__setattr__(hm, "Mmin", Mmin)
            object.__setattr__(hm, "M0"  , Mmin)
            object.__setattr__(hm, "M1"  , M1  )
            galaxyDensity.append( hm.averageGalaxyDensity() )
        galaxyDensities.append( np.array(galaxyDensity) )

    fig = plt.figure()
    plt.loglog()
    for galaxyDensity, M1 in zip(galaxyDensities, M1List):
        plt.plot(MminList, galaxyDensity, "-.", zorder = 2, label = f"{M1:.2e}")
    plt.axhline(2.5e-04, color = "black", zorder = 1)
    plt.legend(title = "M1 [Msun]", ncols = 3)
    plt.xlabel("Mmin [Msun]")
    plt.ylabel("ngal")
    plt.grid()
    fig.savefig("../test-data/images/galaxy_density.png")
    # plt.show()
    return

def testHaloModel():

    files, unitMass, hubble, hm = setupLoading( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )
    object.__setattr__(hm, "alpha", 1.)

    m = np.logspace(8, 16, 21)
    Ncen = hm.centralCount( np.log(m) )
    Nsat = hm.satelliteCount( np.log(m) )

    fig = plt.figure()
    plt.semilogx()
    plt.plot(m, Ncen, "-.", label = "Ncen")
    plt.plot(m, Nsat, "-.", label = "Nsat")
    plt.legend(ncols = 3)
    plt.xlabel("M [Msun]")
    plt.ylabel("N")
    plt.ylim(-0.1, 5.)
    plt.grid()
    fig.savefig("../test-data/images/halo_model.png")
    # plt.show()
    return

if __name__ == "__main__":
    # testGalaxyCatalogGeneration()
    testMinimizeHaloModel()
    # testSatelliteFraction()
    # testGalaxydensity()
    # testHaloModel()