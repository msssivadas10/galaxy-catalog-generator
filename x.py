
import numpy as np
import matplotlib.pyplot as plt
from haloutils import Eisenstein98_zb, HaloModel, halo_massfunction, halo_bias

ps = Eisenstein98_zb(
    z    =  0. , H0 = 70., Om0    =  0.3, Ob0   = 0.05, 
    Ode0 =  0.7, ns =  1., sigma8 =  1. , Tcmb0 = 2.725, 
    w0   = -1. , wa =  0.,
)
hm = HaloModel(
    ps, 
    lnm_min    = np.log(2.563e+12),
    sigma_m    = 0.,
    lnm0       = np.log(2.563e+12),
    lnm1       = np.log(6.147e+13),
    alpha      = 1.  ,
    scale_shmf = 0.5 ,
    slope_shmf = 2.  ,
    Delta_m    = 200.,
)

def test_linear_growth():
    z = np.linspace(0., 10, 11)
    y = ps.linear_growth(z, 0)
    plt.plot(z, y, 'o-')
    plt.show()

def test_power_spectrum():
    k = np.logspace(-4, 3, 11)
    y = ps.power( np.log(k) )
    plt.loglog(k, y, '-o')
    plt.show()

def test_variance():
    r = np.logspace(-3, 2, 11)
    s = ps.spectral_moment( np.log(r), 0, 0, 'tophat' )
    plt.loglog(r, s, '-o')
    plt.show()

def test_correlation():
    r = np.logspace(-3, 2, 11)
    c = ps.correlation( np.log(r), 0 )
    plt.loglog(r, c, '-o')
    plt.show()

def test_mass_function_models():
    m = np.logspace(6, 16, 21)
    f = np.zeros_like(m)
    for model in ['press74', 'sheth01', 'jenkins01', 'reed03', 'tinker08', 'courtin10', 'crocce10']:
        f, s, _ = halo_massfunction( np.log(m), model, ps, 200., retval = 'fs' )
        plt.semilogx(s, f, 'o-', label = model)
    plt.legend(title = 'halo massfunction')
    plt.show()

def test_mass_function_tinker08():
    m = np.logspace(8, 16, 11)
    n, *_ = halo_massfunction( np.log(m), 'tinker08', ps, 200., retval = 'dn/dm' )
    plt.loglog(m, n, 'o-')
    plt.show()

def test_bias_models():
    m = np.logspace(6, 16, 21)
    f = np.zeros_like(m)
    for model in ['cole89', 'sheth01', 'tinker10']:
        f, s, _ = halo_bias( np.log(m), model, ps, 200. )
        plt.loglog(1.686/s, f, 'o-', label = model)
    plt.legend(title = 'halo bias')
    plt.show()

def test_galaxy_count():
    m = np.logspace(8, 14, 21)
    nc = hm.central_count(np.log(m))
    ns = hm.satellite_count(np.log(m))
    plt.semilogx(m, nc, 'o-', m, ns, 'o-')
    plt.show()

def test_shmf():
    x = np.linspace(0., hm.scale_shmf+0.1, 21)
    y = hm.subhalo_massfunction(x, np.log(1e+13))
    plt.plot(x, y, 'o-')
    plt.show()

def test_halo_concentration():
    m = np.logspace(8, 16, 21)
    c = hm.halo_concentration(np.log(m))
    plt.semilogx(m, c, 'o-')
    plt.show()

def test_galaxy_generation():
    lnm   = np.log(1e+14)
    data  = hm.generate_galaxies(lnm, count = 200_000)
    r_gal = np.sqrt(np.sum(np.square(data[1:,0:3]), axis=1))
    x_gal = data[1:,3]

    r_vir = np.exp( hm.lagrangian_r(lnm) )
    c     = hm.halo_concentration(lnm)
    rs    = r_vir / c
    bins  = np.logspace(np.log10(1e-3*r_vir), np.log10(r_vir), 50)
    hist, edges = np.histogram(r_gal, bins=bins, density=True)
    centers = 0.5*(edges[1:] + edges[:-1])
    rho_theory = 1/(centers/rs*(1+centers/rs)**2)
    rho_theory /= np.trapz(rho_theory*centers**2, centers) # normalize
    plt.loglog(centers, hist, label='Sampled')
    plt.loglog(centers, rho_theory*centers**2, label='Theory')
    plt.show()

    m_halo = np.exp(lnm)
    bins   = np.logspace(np.log10(1e-3*m_halo), np.log10(m_halo), 50)
    hist, edges = np.histogram(x_gal, bins=bins, density=True)
    centers = 0.5*(edges[1:] + edges[:-1])
    rho_theory = hm.subhalo_massfunction(centers/m_halo, np.log(m_halo))
    rho_theory /= np.trapz(rho_theory, centers) # normalize
    plt.loglog(centers, hist, label='Sampled')
    plt.loglog(centers, rho_theory, label='Theory')
    plt.show()

def test_galaxy_generation2():
    from haloutils._core import lib, i8, i4, hmargs_t, POINTER, f8pointer

    halodata_t = np.dtype([('id', 'i8'), ('pos', 'f8', 3), ('mass', 'f8')])
    lib.generate_galaxy_catalog.argtypes = [
        np.ctypeslib.ndpointer(dtype=halodata_t, ndim=1, flags='C_CONTIGUOUS'), 
        i8, f8pointer(2), i8, POINTER(hmargs_t), f8pointer(2), i8, i8, i4,      
    ]
    lib.generate_galaxy_catalog.restype = None

    lib.generate_cspline.argtypes = [ i8, f8pointer(2) ]
    lib.generate_cspline.restype  = None
  
    sigma      = np.zeros((101, 3))
    sigma[:,0] = np.linspace(np.log(1e+08), np.log(1e+16), sigma.shape[0])
    sigma[:,2] = hm.lagrangian_r(sigma[:,0])
    sigma[:,1] = ps.spectral_moment(sigma[:,2])
    sigma[:,1] = np.log(sigma[:,1]) / 2.
    lib.generate_cspline(sigma.shape[0], sigma)

    bbox = np.array([[-100., -100., -100.], 
                     [ 100.,  100.,  100.]])

    rng = np.random.default_rng(12345)

    m_halo = 1e+14
    hcat   = np.zeros((20_000,), dtype=halodata_t)
    hcat['id']   = np.arange(hcat.shape[0]) + 1
    hcat['mass'] = m_halo
    hcat['pos']  = rng.uniform(bbox[0,:], bbox[1,:], (hcat.shape[0], 3))
    # print(hcat)

    fid = 0
    lib.generate_galaxy_catalog( hcat, hcat.shape[0], sigma, sigma.shape[0], hm.pointer(), bbox, 123456, fid, 4 )

    galdata_t = np.dtype([('id', '<i8'), ('pos', '<f8', 3), ('mass', '<f8'), ('typ', 'S1')])
    gcat = np.fromfile(f'{fid}.dat', dtype=galdata_t)
    # print(gcat)
    # print(gcat.shape)

    mask      = gcat['typ'] == b'c'
    central   = gcat[mask]
    satellite = gcat[~mask]

    halo_map   = { key: index for index, key in enumerate(central['id']) }
    halo_order = [ halo_map[key] for key in satellite['id'] ]
    bsize  = bbox[1,:] - bbox[0,:]
    rdist  = central[halo_order]['pos'] - satellite['pos']
    rdist -= bsize * np.rint(rdist / bsize)
    rdist  = np.sqrt(np.sum(rdist**2, axis=-1))
    r_vir  = np.exp( hm.lagrangian_r(np.log(m_halo)) )
    c      = hm.halo_concentration(np.log(m_halo))
    rs     = r_vir / c
    bins   = np.logspace(np.log10(1e-3*r_vir), np.log10(r_vir), 50)
    hist, edges = np.histogram(rdist, bins=bins, density=True)
    centers     = 0.5*(edges[1:] + edges[:-1])
    rho_theory  = 1/(centers/rs*(1+centers/rs)**2)
    rho_theory /= np.trapz(rho_theory*centers**2, centers) # normalize
    plt.loglog(centers, hist, label='Sampled')
    plt.loglog(centers, rho_theory*centers**2, label='Theory')
    plt.show()

    msat = satellite['mass'] 
    bins = np.logspace(np.log10(1e-3*m_halo), np.log10(m_halo), 50)
    hist, edges = np.histogram(msat, bins=bins, density=True)
    centers = 0.5*(edges[1:] + edges[:-1])
    rho_theory = hm.subhalo_massfunction(centers/m_halo, np.log(m_halo))
    rho_theory /= np.trapz(rho_theory, centers) # normalize
    plt.loglog(centers, hist, label='Sampled')
    plt.loglog(centers, rho_theory, label='Theory')
    plt.show()

if __name__ == '__main__':
    test_linear_growth()
    test_power_spectrum()
    test_variance()
    test_correlation()
    test_mass_function_models()
    test_mass_function_tinker08()
    test_bias_models()
    test_galaxy_count()
    test_shmf()
    test_halo_concentration()
    test_galaxy_generation()
    test_galaxy_generation2()
