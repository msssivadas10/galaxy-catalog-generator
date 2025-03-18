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

def setupLoading(path):

    files = glob.glob( path )

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

    return files, unitMass, hubble, hm

def testGalaxyCounts():

    files, unitMass, hubble, hm = setupLoading( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )

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

def _testSatellieMassFunction(massH):

    files, unitMass, hubble, hm = setupLoading( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )

    binsize = 0.01
    # massH   = 1e+13

    fracBinEdges = np.linspace(hm.Mmin*hubble / massH, hm.scaleSHMF, 21)
    fracMass     = ( fracBinEdges[1:] + fracBinEdges[:-1] ) / 2.

    satelliteMassFunction = np.zeros(fracMass.shape)
    for file in files:
        with asdf.open(file) as af:

            uniqueIDs, centralIndex, inverseIndex = np.unique( 
                af["data"]["parentHaloID"], 
                return_index   = True, 
                return_inverse = True, 
            )

            # Satellite galaxies
            satelliteIndex, = np.where( af["data"]["galaxyType"] == 2 )

            # Mass of parent halos of satellites 
            haloMass = af["data"]["galaxyMass"][centralIndex] 
            haloMass = haloMass[inverseIndex[satelliteIndex]]
            
            # Mass of the satllites as a fraction of the halo / central galaxy mass
            massFraction = af["data"]["galaxyMass"][satelliteIndex] / haloMass
            
            # Apply mass bin
            minMass = np.exp(np.log(massH) * (1 - binsize / 2))
            maxMass = np.exp(np.log(massH) * (1 + binsize / 2))
            galaxiesInsideBin, = np.where( ( haloMass > minMass ) & ( haloMass < maxMass ) )
            haloMass       = haloMass[galaxiesInsideBin]
            massFraction   = massFraction[galaxiesInsideBin]
            
            _satelliteMassFunc, _     = np.histogram( massFraction, bins = fracBinEdges )
            satelliteMassFunction[:] += _satelliteMassFunc[:]

    satelliteMassFunction /= np.sum( satelliteMassFunction ) * np.diff(fracBinEdges) 

    L, H, a = hm.Mmin*hubble / massH, hm.scaleSHMF, hm.slopeSHMF
    satelliteMassFunctionExp = a*L**a * fracMass**(-a-1) / ( 1 - (L/H)**a )

    return hm, fracMass, satelliteMassFunction, satelliteMassFunctionExp

def testSatellieMassFunction():
    
    fig = plt.figure()
    plt.loglog()

    for massH, colour in  [(5e+12, "green"), (1e+13, "tab:blue"), (5e+13, "red")]:
        hm, fracMass1, satelliteMassFunction1, satelliteMassFunctionExp1 = _testSatellieMassFunction(massH)
        plt.plot( fracMass1, satelliteMassFunction1   , 's' , color = colour )
        plt.plot( fracMass1, satelliteMassFunctionExp1, '-.', color = colour )
        plt.plot( [], [], '-s', color = colour, label = f"{massH:.3g}" )

    plt.legend(ncols = 3)
    plt.xlabel("Galaxy - Halo mass ratio")
    plt.ylabel("PDF")
    
    plt.grid()
    # plt.show()
    fig.savefig("../test-data/images/galaxy_mass_pdf_hugebase_z3.png")

    return

def _testHaloProfile(massH):

    files, unitMass, hubble, hm = setupLoading( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )

    binsize = 0.01
    # massH   = 1e+13
    rho_h   = hm.mfmodel.rho_m * (1 + hm.redshift)**3
    radiusH = np.cbrt( 0.75 / np.pi * ( (massH / hubble) / rho_h ) ) * hubble # halo radius in Mpc/h
    concH   = hm.haloConcentration( np.log(massH / hubble) )
    Rs      = radiusH / concH
    # print( concH )

    distBinEdges = np.logspace(-3, np.log10(radiusH), 21)
    distGH       = np.sqrt( distBinEdges[1:] * distBinEdges[:-1] )

    haloProfile = np.zeros(distGH.shape)
    for file in files:
        with asdf.open(file) as af:

            uniqueIDs, centralIndex, inverseIndex = np.unique( 
                af["data"]["parentHaloID"], 
                return_index   = True, 
                return_inverse = True, 
            )

            # Mass of parent halos of galaxies 
            haloMass = af["data"]["galaxyMass"][centralIndex][inverseIndex]

            # Distance of the galaxy from the halo center
            posX, posY, posZ = np.transpose( af["data"]["galaxyPosition"] )
            deltaX           = posX - posX[centralIndex][inverseIndex]
            deltaY           = posY - posY[centralIndex][inverseIndex]
            deltaZ           = posZ - posZ[centralIndex][inverseIndex]
            galaxyDistance   = np.sqrt( deltaX**2 + deltaY**2 + deltaZ**2 )

            # Apply mass bin
            minMass = np.exp(np.log(massH) * (1 - binsize / 2))
            maxMass = np.exp(np.log(massH) * (1 + binsize / 2))
            galaxiesInsideBin, = np.where( ( haloMass > minMass ) & ( haloMass < maxMass ) )
            haloMass       = haloMass[galaxiesInsideBin]
            galaxyDistance = galaxyDistance[galaxiesInsideBin]
            
            _haloProfile, _  = np.histogram( galaxyDistance, bins = distBinEdges )
            haloProfile[:]  += _haloProfile[:]

    haloProfile /= np.sum( haloProfile ) * np.diff(distBinEdges) 

    x     = distGH / Rs
    haloProfileExp = x / (x + 1)**2 / ( np.log(1 + concH) - concH / ( 1 + concH ) ) / Rs

    return hm, distGH, haloProfile, haloProfileExp

def testHaloProfile():

    fig = plt.figure()
    plt.loglog()

    for massH, colour in  [(5e+12, "green"), (1e+13, "tab:blue"), (5e+13, "red")]:
        hm, distGH, haloProfile, haloProfileExp = _testHaloProfile(massH)
        plt.plot( distGH, haloProfile   , 's' , color = colour )
        plt.plot( distGH, haloProfileExp, '-.', color = colour )
        plt.plot( [], [], '-s', color = colour, label = f"{massH:.3g}" )

    plt.legend(ncols = 3)
    plt.xlabel("distance [Mpc/h]")
    plt.ylabel("$r^2 \\rho(r)$")
    
    plt.grid()
    # plt.show()
    fig.savefig("../test-data/images/halo_profile_nfw_hugebase_z3.png")

    return

if __name__ == "__main__":
    # testGalaxyCatalogGeneration()
    # testGalaxyCounts()
    # testSatellieMassFunction()
    testHaloProfile()

# a, b = 0.1, 0.5
# p = 1.91
# u = np.random.uniform(size = 10000)
# x = ( -( u*b**p - u*a**p - b**p ) / ( a*b )**p )**( -1./p )
# _, xb, _ = plt.hist(x, bins=21, density=True)
# xb = np.linspace(a, b, 51)
# yb = p * a**p * xb**(-p-1) / (1 - (a/b)**p)
# plt.plot(xb, yb, '--')
# plt.show()

# c = np.logspace(-3, 1.5, 21)
# A = np.log(1 + c) - c / (1 + c)
# # plt.loglog(c, A)
# # plt.show()

# from scipy.interpolate import CubicSpline
# ctable = CubicSpline(A, c)

# c = 5.31
# u = np.random.uniform(size = 10000)
# x = ctable( u * ( np.log(1 + c) - c / (1 + c) ) )
# _, xb, _ = plt.hist(x, bins=21, density=True)
# xb = ( xb[1:] + xb[:-1] ) / 2
# yb = xb / (1 + xb)**2 / ( np.log(1 + c) - c / (1 + c) )
# plt.plot(xb, yb, '--')
# plt.show()

# files, unitMass, hubble, hm = setupLoading( "../test-data/galaxy_catalog/z3.000/galaxy_info_*.asdf" )

# massH = 1e+13
# rho_h   = hm.mfmodel.rho_m * (1 + hm.redshift)**3
# radiusH = np.cbrt( 0.75 / np.pi * ( (massH / hubble) / rho_h ) ) * hubble # halo radius in Mpc/h
# concH   = hm.haloConcentration( np.log(massH / hubble) )
# Rs      = radiusH / concH

# posH = np.random.uniform(-500, 500, (10000, 3))
# satellitePos = []
# for id, _posH in enumerate(posH):
#     posG = hm.generateSatellitePositions(np.log(massH / hubble), _posH / hubble) 
#     if posG is not None:
#         satellitePos.append( np.hstack([ [[id]] * posG.shape[0], posG * hubble ]) )
# satellitePos = np.vstack( satellitePos )

# uniqueIDs, centralIndex, inverseIndex = np.unique( 
#     satellitePos[:,0], 
#     return_index   = True, 
#     return_inverse = True, 
# )

# x = np.sqrt( np.sum( (satellitePos[:,1:4] - satellitePos[centralIndex[inverseIndex],1:4])**2, axis = 1 ) )
# # print( max(x), radiusH )
# # x = satellitePos

# plt.loglog()
# y, xb = np.histogram(x, bins=np.logspace(-3, np.log10(radiusH), 21), density=False)
# y  = y / np.sum( y ) / ( np.diff(xb) )
# xb = ( xb[1:] + xb[:-1] ) / 2
# yb = (xb/Rs) / (1 + (xb/Rs))**2 / ( np.log(1 + concH) - concH / (1 + concH) ) / Rs
# plt.plot(xb, y , 's')
# plt.plot(xb, yb, '-.')
# plt.show()



