
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
args_ps.Obh2    = 0.04        * args_ps.h**2    
args_ps.ns      = 0.9  
args_ps.sigma8  = 0.8     
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

pktab = np.zeros((1001, 2))
pktab[:, 0] = np.linspace( np.log(1e-04), np.log(1e+04), pktab.shape[0] )
pktab[:, 1] = ps_eisenstein98_zb(pktab[:, 0], ctypes.byref(args_ps))

var8_calc  = variance( np.log(8./args_ps.h), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] )
args_ps.norm  = args_ps.sigma8**2 / var8_calc 
pktab[:,1] = pktab[:,1] + np.log(args_ps.norm)
# print( variance( np.log(8./args_ps.h), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] ), args_ps.sigma8**2 )

args_mf = hmfargs_t()
args_mf.z       = args_ps.z
args_mf.lnm     = np.log(1e+10)
args_mf.H0      = 70.0
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

# z = np.linspace(0., 10, 11)
# y = linear_growth(z, 0, ctypes.byref(args_gf))
# plt.plot(z, y, 'o-')
# plt.show()

# k = np.logspace(-4, 3, 11)
# y = np.exp( ps_eisenstein98_zb(np.log(k), ctypes.byref(args_ps)) )
# plt.loglog(k, y, '-o')
# plt.show()

# r = np.logspace(-3, 2, 11)
# s = variance( np.log(r), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] )
# c = correlation( np.log(r), 0, pktab, pktab.shape[0], pktab.shape[1] )
# plt.loglog(r, s, '-o')
# plt.loglog(r, c, '-o')
# plt.show()

# s = np.logspace(-1, 3, 21)
# f = np.zeros_like(s)
# for model in ['press74', 'sheth01', 'jenkins01', 'reed03', 'tinker08', 'courtin10', 'crocce10']:
#     func = getattr(lib, 'hmf_' + model)
#     for i in range(s.shape[0]):
#         args_mf.s = s[i]
#         f[i]      = func(ctypes.byref(args_mf), -1) 
#     plt.semilogx(s, f, 'o-', label = model)
# plt.legend(title = 'halo massfunction')
# plt.show()

# s = np.logspace(-1, 3, 21)
# f = np.zeros_like(s)
# for model in ['cole89', 'sheth01', 'tinker10']:
#     func = getattr(lib, 'hbf_' + model)
#     for i in range(s.shape[0]):
#         args_mf.s = s[i]
#         f[i]      = func(ctypes.byref(args_mf)) 
#     plt.loglog(1.686/s, f, 'o-', label = model)
# plt.legend(title = 'halo bias')
# plt.show()

# m = np.logspace(8, 20, 11)
# n = np.zeros_like(m)
# for i in range(m.shape[0]):
#     args_mf.lnm = np.log( m[i] )
#     lib.setup_hmf_tinker08(ctypes.byref(args_mf), 0, pktab, pktab.shape[0], pktab.shape[1])
#     n[i] = lib.hmf_tinker08(ctypes.byref(args_mf), 0)
# plt.loglog(m, n, 'o-')
# plt.show()