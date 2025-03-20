import warnings
warnings.catch_warnings(action="ignore")

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from galaxy_catalog_generator.halomodel import HaloModel

hm = HaloModel.create(
    Mmin      = 1e+12, 
    sigmaM    = 0., 
    M0        = 1e+12, 
    M1        = 1e+12, 
    alpha     = 0.3, 
    scaleSHMF = 0.50, 
    slopeSHMF = 2.00, 
    redshift  = 0.,
    cosmo     = FlatLambdaCDM(H0 = 70., Om0 = 0.3, Ob0 = 0.05, Tcmb0 = 2.725), 
    psmodel   = "eisenstein98_zb", 
    mfmodel   = "tinker08", 
    ns        = 1., 
    sigma8    = 1., 
    Delta     = 200,
)

def testGalaxyCatalogGeneration():
    # Test galaxy catalog generation

    haloMass = 1e+13 # in Msun
    satelliteCatalog = hm.generateSatellitePositions( np.log(haloMass), [10., 10., 10.] )
    print(satelliteCatalog)
    return

# m = np.logspace(7, 16, 201)
# y = hm.haloMassFunction( np.log(m), "dndlnm" )
# y1 = hm.centralCount( np.log(m) )   * y
# y2 = hm.satelliteCount( np.log(m) ) * y

# plt.loglog()
# plt.plot(m, y , "-.", color = "tab:blue")
# plt.plot(m, y1, "-.", color = "red")
# plt.plot(m, y2, "-.", color = "orange")
# plt.show()

# print( hm.averageHaloDensity() )
# print( hm.averageGalaxyDensity() )
# print( hm.averageSatelliteFraction() )

# x = np.logspace(10, 18, 51)
# y = 10**( 3.16 * np.log10(x) - 24.3 )

# plt.loglog()
# plt.plot(x, y , "-.", color = "tab:blue")
# plt.show()

observedGalaxyDensity = 2.5e-04
observedSatelliteFraction = 0.1

# from scipy.optimize import minimize

# def scoreFunction(m):
#     object.__setattr__(hm, "Mmin", m[0])
#     object.__setattr__(hm, "M0"  , m[0])
#     object.__setattr__(hm, "M1"  , m[1]  )
#     galaxyDensity     = hm.averageGalaxyDensity()
#     satelliteFraction = hm.averageSatelliteFraction()
#     score = (galaxyDensity - observedGalaxyDensity) **2 + (satelliteFraction - observedSatelliteFraction)**2
#     return score

# res = minimize(
#     scoreFunction, 
#     x0 = [1e+11, 1e+11], 
#     bounds = [(1e+10, 1e+18), (1e+10, 1e+18)]
# )
# print( res )

x = np.logspace(10, 17, 51)
Mmin, M1 = np.meshgrid(x, x)
score = np.zeros_like(Mmin)
for i in range( score.shape[0] ):
    for j in range( score.shape[1] ):
        print(i, j)
        object.__setattr__(hm, "Mmin", Mmin[i, j])
        object.__setattr__(hm, "M0"  , Mmin[i, j])
        object.__setattr__(hm, "M1"  , M1[i, j]  )
        galaxyDensity     = hm.averageGalaxyDensity()
        satelliteFraction = hm.averageSatelliteFraction()
        score[i, j] = (galaxyDensity - observedGalaxyDensity) **2 + (satelliteFraction - observedSatelliteFraction)**2

i, j = np.unravel_index(np.argmin(score), score.shape)
optMmin = Mmin[i, j]
optM1   = M1[i, j]

plt.loglog()
plt.contourf(Mmin, M1, np.log10(score), levels = 21)
plt.plot([optMmin], [optM1], 'o', color = "black")
plt.text(optMmin, optM1, f"  ({optMmin:2e}, {optM1:.2e})")
plt.colorbar()
plt.xlabel("Mmin [Msun]")
plt.ylabel("M1 [Msun]")
plt.show()
