
import numpy as np
import matplotlib.pyplot as plt
import numpy.ctypeslib as npct
import ctypes

# load the shared library
lib = npct.load_library("libpowerspectrum", ".")

class zfargs_t(ctypes.Structure):
    _fields_ = [
        ('Om0'    , ctypes.c_double), 
        ('Ode0'   , ctypes.c_double),
        ('w0'     , ctypes.c_double),
        ('wa'     , ctypes.c_double),
        ('abstol' , ctypes.c_double),
        ('reltol' , ctypes.c_double),
        ('maxiter', ctypes.c_int64 ),
    ]

lib.linear_growth.argtypes = [
    ctypes.c_double, 
    ctypes.c_int, 
    ctypes.POINTER(zfargs_t), 
]
lib.linear_growth.restype = ctypes.c_double
linear_growth = np.vectorize(lib.linear_growth, excluded=[2])

class psargs_t(ctypes.Structure):
    _fields_ = [
        ('z'         , ctypes.c_double  ),
        ('h'         , ctypes.c_double  ), 
        ('Omh2'      , ctypes.c_double  ),
        ('Obh2'      , ctypes.c_double  ),
        ('Onuh2'     , ctypes.c_double  ),
        ('Nnu'       , ctypes.c_double  ),        
        ('ns'        , ctypes.c_double  ),
        ('sigma8'    , ctypes.c_double  ),
        ('theta'     , ctypes.c_double  ),
        ('dplus_z'   , ctypes.c_double  ),        
        ('dplus_0'   , ctypes.c_double  ),        
        ('norm'      , ctypes.c_double  ),        
        ('include_nu', ctypes.c_int     ),
        ('z_eq'      , ctypes.c_double  ),
        ('z_d'       , ctypes.c_double  ),        
        ('s'         , ctypes.c_double  ),
        ('k_silk'    , ctypes.c_double  ),
        ('param'     , ctypes.c_double*5),
    ]

CALLBACK = ctypes.CFUNCTYPE(
    ctypes.c_double,    # return type
    ctypes.c_double,    # x input
    ctypes.c_void_p     # ctx pointer
)

def as_callback(fn):
    def wrapper(lnk, argsp):
        args_ps = ctypes.cast(argsp, ctypes.POINTER(psargs_t)).contents 
        return fn(lnk, ctypes.byref(args_ps))
    return CALLBACK(wrapper)

# define argument and return types
lib.init_eisenstein98_zb.argtypes = [ ctypes.POINTER(psargs_t) ]
lib.init_eisenstein98_zb.restype  = None

lib.ps_eisenstein98_zb.argtypes = [ ctypes.c_double, ctypes.POINTER(psargs_t) ]
lib.ps_eisenstein98_zb.restype  = ctypes.c_double
ps_eisenstein98_zb = np.vectorize(lib.ps_eisenstein98_zb)

lib.variance.argtypes = [
    ctypes.c_double, 
    ctypes.c_int, 
    ctypes.c_int, 
    ctypes.c_int, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64,
    ctypes.c_int,
]
lib.variance.restype = ctypes.c_double
variance = np.vectorize(lib.variance, excluded=[4])

lib.correlation.argtypes = [
    ctypes.c_double, 
    ctypes.c_int, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64,
    ctypes.c_int,
]
lib.correlation.restype = ctypes.c_double
correlation = np.vectorize(lib.correlation, excluded=[2])

class hmfargs_t(ctypes.Structure):
    _fields_ = [
        ('z'       , ctypes.c_double   ),
        ('lnm'     , ctypes.c_double   ),
        ('H0'      , ctypes.c_double   ),        
        ('Om0'     , ctypes.c_double   ),
        ('Delta_m' , ctypes.c_double   ),
        ('s'       , ctypes.c_double   ),
        ('dlnsdlnm', ctypes.c_double   ), 
        ('rho_m'   , ctypes.c_double   ),
        ('param'   , ctypes.c_double*16),        
    ]

lib.setup_hmf.argtypes = [
    ctypes.POINTER(hmfargs_t), 
    ctypes.c_int, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64,
    ctypes.c_int,
]
lib.setup_hmf.restype = None

lib.setup_hmf_tinker08.argtypes = [
    ctypes.POINTER(hmfargs_t), 
    ctypes.c_int, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64,
    ctypes.c_int,
]
lib.setup_hmf_tinker08.restype = None

lib.hmf_press74.argtypes = [ ctypes.POINTER(hmfargs_t), ctypes.c_int ]
lib.hmf_press74.restype  = ctypes.c_double

lib.hmf_sheth01.argtypes = [ ctypes.POINTER(hmfargs_t), ctypes.c_int ]
lib.hmf_sheth01.restype  = ctypes.c_double

lib.hmf_jenkins01.argtypes = [ ctypes.POINTER(hmfargs_t), ctypes.c_int ]
lib.hmf_jenkins01.restype  = ctypes.c_double

lib.hmf_reed03.argtypes = [ ctypes.POINTER(hmfargs_t), ctypes.c_int ]
lib.hmf_reed03.restype  = ctypes.c_double

lib.hmf_tinker08.argtypes = [ ctypes.POINTER(hmfargs_t), ctypes.c_int ]
lib.hmf_tinker08.restype  = ctypes.c_double

lib.hmf_courtin10.argtypes = [ ctypes.POINTER(hmfargs_t), ctypes.c_int ]
lib.hmf_courtin10.restype  = ctypes.c_double

lib.hmf_crocce10.argtypes = [ ctypes.POINTER(hmfargs_t), ctypes.c_int ]
lib.hmf_crocce10.restype  = ctypes.c_double

lib.hbf_cole89.argtypes = [ ctypes.POINTER(hmfargs_t) ] 
lib.hbf_cole89.restype  = ctypes.c_double

lib.hbf_sheth01.argtypes = [ ctypes.POINTER(hmfargs_t) ] 
lib.hbf_sheth01.restype  = ctypes.c_double

lib.hbf_tinker10.argtypes = [ ctypes.POINTER(hmfargs_t) ]
lib.hbf_tinker10.restype  = ctypes.c_double

class hmargs_t(ctypes.Structure):
    _fields_ = [
        ('lnm_min'   , ctypes.c_double),
        ('sigma_m'   , ctypes.c_double),
        ('lnm0'      , ctypes.c_double),
        ('lnm1'      , ctypes.c_double),
        ('alpha'     , ctypes.c_double),
        ('scale_shmf', ctypes.c_double),
        ('slope_shmf', ctypes.c_double),
        ('z'         , ctypes.c_double),
        ('H0'        , ctypes.c_double),
        ('Om0'       , ctypes.c_double),
        ('Om0'       , ctypes.c_double), 
        ('Ode0'      , ctypes.c_double),
        ('w0'        , ctypes.c_double),
        ('wa'        , ctypes.c_double),
        ('Delta_m'   , ctypes.c_double),
        ('dplus'     , ctypes.c_double),
    ]        

class cgargs_t(ctypes.Structure):
    _fields_ = [
        ('lnm'    , ctypes.c_double  ), 
        ('pos'    , ctypes.c_double*3),
        ('lnr'    , ctypes.c_double  ),
        ('s'      , ctypes.c_double  ),
        ('c'      , ctypes.c_double  ),
        ('n_cen'  , ctypes.c_int64   ),
        ('n_sat'  , ctypes.c_int64   ),
        ('rstate' , ctypes.c_int64   ),
        ('boxsize', ctypes.c_double*3),
        ('offset' , ctypes.c_double*3),
    ]

lib.lagrangian_r.argtypes = [ ctypes.POINTER(hmargs_t), ctypes.c_double ]
lib.lagrangian_r.restype  = ctypes.c_double
lagrangian_r = np.vectorize(lib.lagrangian_r, excluded=[0])

lib.central_count.argtypes = [ ctypes.POINTER(hmargs_t), ctypes.c_double ]
lib.central_count.restype  = ctypes.c_double
central_count = np.vectorize(lib.central_count, excluded=[0])

lib.satellite_count.argtypes = [ ctypes.POINTER(hmargs_t), ctypes.c_double ]
lib.satellite_count.restype  = ctypes.c_double
satellite_count = np.vectorize(lib.satellite_count, excluded=[0])

lib.subhalo_mass_function.argtypes = [ ctypes.POINTER(hmargs_t), ctypes.c_double, ctypes.c_double ]
lib.subhalo_mass_function.restype  = ctypes.c_double
subhalo_mass_function = np.vectorize(lib.subhalo_mass_function, excluded=[0, 2])

lib.halo_concentration.argtypes = [ ctypes.POINTER(hmargs_t), ctypes.c_double ]
lib.halo_concentration.restype  = ctypes.c_double
halo_concentration = np.vectorize(lib.halo_concentration, excluded=[0])

lib.setup_catalog_generation.argtypes = [
    ctypes.POINTER(hmargs_t),
    ctypes.POINTER(cgargs_t),
]
lib.setup_catalog_generation.restype  = None

lib.generate_galaxies.argtypes = [
    ctypes.POINTER(hmargs_t),
    ctypes.POINTER(cgargs_t),
    ctypes.c_int64,
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
]
lib.generate_galaxies.restype = None

lib.average_halo_density.argtypes = [
    ctypes.POINTER(hmargs_t), 
    ctypes.c_double,         
    ctypes.c_double, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64, 
    ctypes.c_double, 
    ctypes.c_double, 
    ctypes.c_int64,
]
lib.average_halo_density.restype = ctypes.c_double

lib.average_galaxy_density.argtypes = [
    ctypes.POINTER(hmargs_t), 
    ctypes.c_double,         
    ctypes.c_double, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64, 
    ctypes.c_double, 
    ctypes.c_double, 
    ctypes.c_int64,
]
lib.average_galaxy_density.restype = ctypes.c_double

lib.average_satellite_frac.argtypes = [
    ctypes.POINTER(hmargs_t), 
    ctypes.c_double,         
    ctypes.c_double, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64, 
    ctypes.c_double, 
    ctypes.c_double, 
    ctypes.c_int64,
]
lib.average_satellite_frac.restype = ctypes.c_double

lib.average_galaxy_bias.argtypes = [
    ctypes.POINTER(hmargs_t), 
    ctypes.c_double,         
    ctypes.c_double, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64, 
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), 
    ctypes.c_int64, 
    ctypes.c_double, 
    ctypes.c_double, 
    ctypes.c_int64,
]
lib.average_galaxy_bias.restype = ctypes.c_double

############################################################################################################

args_gf = zfargs_t()
args_gf.Om0     =  0.3
args_gf.Ode0    =  0.7
args_gf.w0      = -1.
args_gf.wa      =  0.
args_gf.abstol  =  1e-08
args_gf.reltol  =  1e-08
args_gf.maxiter =  100

args_ps = psargs_t()
args_ps.z       = 0. 
args_ps.h       = 0.7 
args_ps.Omh2    = args_gf.Om0 * args_ps.h**2    
args_ps.Obh2    = 0.05        * args_ps.h**2    
args_ps.ns      = 1.0  
args_ps.sigma8  = 1.0     
args_ps.theta   = 2.725 / 2.7    
args_ps.norm    = 1.    
args_ps.dplus_z = linear_growth(args_ps.z, 0, ctypes.byref(args_gf))       
args_ps.dplus_0 = linear_growth(0.       , 0, ctypes.byref(args_gf)) 

lib.init_eisenstein98_zb( ctypes.byref(args_ps) )
# print('z          = ', args_ps.z         )
# print('h          = ', args_ps.h         )
# print('Omh2       = ', args_ps.Omh2      )
# print('Obh2       = ', args_ps.Obh2      )
# print('Onuh2      = ', args_ps.Onuh2     )
# print('Nnu        = ', args_ps.Nnu       )
# print('ns         = ', args_ps.ns        )
# print('sigma8     = ', args_ps.sigma8    )
# print('theta      = ', args_ps.theta     )
# print('dplus_z    = ', args_ps.dplus_z   )
# print('dplus_0    = ', args_ps.dplus_0   )
# print('norm       = ', args_ps.norm      )
# print('include_nu = ', args_ps.include_nu)
# print('z_eq       = ', args_ps.z_eq      )
# print('z_d        = ', args_ps.z_d       )
# print('s          = ', args_ps.s         )
# print('k_silk     = ', args_ps.k_silk    )
# print('param      = ', list(args_ps.param))

pktab        = np.zeros((1001, 2))
pktab[:, 0]  = np.linspace( np.log(1e-04), np.log(1e+04), pktab.shape[0] )
pktab[:, 1]  = ps_eisenstein98_zb(pktab[:, 0], ctypes.byref(args_ps))
var8_calc    = variance( np.log(8./args_ps.h), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] )
args_ps.norm = args_ps.sigma8**2 / var8_calc 
pktab[:,1]   = pktab[:,1] + np.log(args_ps.norm)

args_mf = hmfargs_t()
args_mf.z       = args_ps.z
args_mf.lnm     = np.log(1e+10)
args_mf.H0      = 100.0*args_ps.h
args_mf.Om0     = args_gf.Om0
args_mf.Delta_m = 200.

lib.setup_hmf_tinker08(ctypes.byref(args_mf), 0, pktab, pktab.shape[0], pktab.shape[1])
# print('z        = ', args_mf.z          )
# print('lnm      = ', args_mf.lnm        )
# print('H0       = ', args_mf.H0         )
# print('Om0      = ', args_mf.Om0        )
# print('Delta_m  = ', args_mf.Delta_m    )
# print('s        = ', args_mf.s          )
# print('dlnsdlnm = ', args_mf.dlnsdlnm   )
# print('rho_m    = ', args_mf.rho_m      )
# print('param    = ', list(args_mf.param))

args_hm = hmargs_t()
args_hm.lnm_min    = np.log(2.563e+12)
args_hm.sigma_m    = 0.
args_hm.lnm0       = np.log(2.563e+12)
args_hm.lnm1       = np.log(6.147e+13)
args_hm.alpha      = 1.
args_hm.scale_shmf = 0.5
args_hm.slope_shmf = 2.
args_hm.z          = args_ps.z
args_hm.H0         = 100.0*args_ps.h
args_hm.Om0        = args_gf.Om0
args_hm.Ode0       = args_gf.Ode0
args_hm.w0         = args_gf.w0  
args_hm.wa         = args_gf.wa  
args_hm.Delta_m    = 200.
args_hm.dplus      = (args_ps.dplus_z / args_ps.dplus_0)

def test_linear_growth():
    z = np.linspace(0., 10, 11)
    y = linear_growth(z, 0, ctypes.byref(args_gf))
    plt.plot(z, y, 'o-')
    plt.show()

def test_power_spectrum():
    k = np.logspace(-4, 3, 11)
    y = np.exp( ps_eisenstein98_zb(np.log(k), ctypes.byref(args_ps)) )
    plt.loglog(k, y, '-o')
    plt.show()

def test_variance():
    r = np.logspace(-3, 2, 11)
    s = variance( np.log(r), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] )
    plt.loglog(r, s, '-o')
    plt.show()

def test_correlation():
    r = np.logspace(-3, 2, 11)
    c = correlation( np.log(r), 0, pktab, pktab.shape[0], pktab.shape[1] )
    plt.loglog(r, c, '-o')
    plt.show()

def test_mass_function_models():
    s = np.logspace(-1, 3, 21)
    f = np.zeros_like(s)
    for model in ['press74', 'sheth01', 'jenkins01', 'reed03', 'tinker08', 'courtin10', 'crocce10']:
        func = getattr(lib, 'hmf_' + model)
        for i in range(s.shape[0]):
            args_mf.s = s[i]
            f[i]      = func(ctypes.byref(args_mf), -1) 
        plt.semilogx(s, f, 'o-', label = model)
    plt.legend(title = 'halo massfunction')
    plt.show()

def test_mass_function_tinker08():
    m = np.logspace(8, 20, 11)
    n = np.zeros_like(m)
    for i in range(m.shape[0]):
        args_mf.lnm = np.log( m[i] )
        lib.setup_hmf_tinker08(ctypes.byref(args_mf), 0, pktab, pktab.shape[0], pktab.shape[1])
        n[i] = lib.hmf_tinker08(ctypes.byref(args_mf), 0)
    plt.loglog(m, n, 'o-')
    plt.show()

def test_bias_models():
    s = np.logspace(-1, 3, 21)
    f = np.zeros_like(s)
    for model in ['cole89', 'sheth01', 'tinker10']:
        func = getattr(lib, 'hbf_' + model)
        for i in range(s.shape[0]):
            args_mf.s = s[i]
            f[i]      = func(ctypes.byref(args_mf)) 
        plt.loglog(1.686/s, f, 'o-', label = model)
    plt.legend(title = 'halo bias')
    plt.show()

def test_galaxy_count():
    m = np.logspace(8, 14, 21)
    nc = central_count(ctypes.byref(args_hm), np.log(m))
    ns = satellite_count(ctypes.byref(args_hm), np.log(m))
    plt.semilogx(m, nc, 'o-', m, ns, 'o-')
    plt.show()

def test_shmf():
    x = np.linspace(0., args_hm.scale_shmf+0.1, 21)
    y = subhalo_mass_function(ctypes.byref(args_hm), x, np.log(1e+13))
    plt.plot(x, y, 'o-')
    plt.show()

def test_halo_concentration():
    m = np.logspace(8, 16, 21)
    r = np.exp( lagrangian_r(ctypes.byref(args_hm), np.log(m)) )
    s = np.sqrt( variance( np.log(r), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] ) )
    c = halo_concentration(ctypes.byref(args_hm), s)
    plt.semilogx(m, c, 'o-')
    plt.show()

def test_galaxy_generation():
    # import time
    # seed = int( time.clock_gettime_ns(0) )
    args_cg = cgargs_t()
    args_cg.lnm = np.log(1e+14)
    args_cg.pos = (0., 0., 0.)
    args_cg.s   = np.sqrt( variance( lagrangian_r( ctypes.byref(args_hm), args_cg.lnm), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] ) )
    args_cg.rstate = 1234567 # seed

    # lib.setup_catalog_generation(ctypes.byref(args_hm), ctypes.byref(args_cg))
    # print('lnm    = ', args_cg.lnm      )
    # print('pos    = ', list(args_cg.pos))
    # print('lnr    = ', args_cg.lnr      )
    # print('s      = ', args_cg.s        )
    # print('c      = ', args_cg.c        )
    # print('n_cen  = ', args_cg.n_cen    )
    # print('n_sat  = ', args_cg.n_sat    )
    # print('rstate = ', args_cg.rstate   )

    n_sat  = 200_000
    r_gal, x_gal = [], [] 
    while n_sat > 0:
        lib.setup_catalog_generation(ctypes.byref(args_hm), ctypes.byref(args_cg))
        n    = args_cg.n_cen+args_cg.n_sat
        data = np.empty((n, 4))
        lib.generate_galaxies(ctypes.byref(args_hm), ctypes.byref(args_cg), n, data)
        r_gal.append( np.sqrt(np.sum(np.square(data[1:,0:3]), axis=1)) )
        x_gal.append( data[1:,3] )
        n_sat -= n
    r_gal = np.hstack(r_gal)
    x_gal = np.hstack(x_gal)

    r_vir = np.exp( args_cg.lnr )
    c     = args_cg.c
    rs    = r_vir / c
    bins  = np.logspace(np.log10(1e-3*r_vir), np.log10(r_vir), 50)
    hist, edges = np.histogram(r_gal, bins=bins, density=True)
    centers = 0.5*(edges[1:] + edges[:-1])
    rho_theory = 1/(centers/rs*(1+centers/rs)**2)
    rho_theory /= np.trapz(rho_theory*centers**2, centers) # normalize
    plt.loglog(centers, hist, label='Sampled')
    plt.loglog(centers, rho_theory*centers**2, label='Theory')
    plt.show()

    m_halo = np.exp(args_cg.lnm)
    bins   = np.logspace(np.log10(1e-3*m_halo), np.log10(m_halo), 50)
    hist, edges = np.histogram(x_gal, bins=bins, density=True)
    centers = 0.5*(edges[1:] + edges[:-1])
    rho_theory = subhalo_mass_function(ctypes.byref(args_hm), centers/m_halo, np.log(m_halo))
    rho_theory /= np.trapz(rho_theory, centers) # normalize
    plt.loglog(centers, hist, label='Sampled')
    plt.loglog(centers, rho_theory, label='Theory')
    plt.show()

if __name__ == '__main__':
    test_linear_growth()
    test_power_spectrum()
    test_correlation()
    test_mass_function_models()
    test_mass_function_tinker08()
    test_bias_models()
    test_galaxy_count()
    test_shmf()
    test_halo_concentration()
    test_galaxy_generation()