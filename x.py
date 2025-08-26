
import numpy as np
import matplotlib.pyplot as plt
import numpy.ctypeslib as npct
import ctypes

# load the shared library
lib = npct.load_library("libpowerspectrum", ".")

class lgargs_t(ctypes.Structure):
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
    ctypes.POINTER(lgargs_t), 
]
lib.linear_growth.restype = ctypes.c_double
linear_growth = np.vectorize(lib.linear_growth, excluded=[2])

args_gf = lgargs_t()
args_gf.Om0     =  0.3
args_gf.Ode0    =  0.7
args_gf.w0      = -1.
args_gf.wa      =  0.
args_gf.abstol  =  1e-08
args_gf.reltol  =  1e-08
args_gf.maxiter =  50

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
        args = ctypes.cast(argsp, ctypes.POINTER(psargs_t)).contents 
        return fn(lnk, ctypes.byref(args))
    return CALLBACK(wrapper)

# define argument and return types
lib.init_eisenstein98_zb.argtypes = [ ctypes.POINTER(psargs_t) ]
lib.init_eisenstein98_zb.restype  = None

lib.ps_eisenstein98_zb.argtypes = [ ctypes.c_double, ctypes.POINTER(psargs_t) ]
lib.ps_eisenstein98_zb.restype  = ctypes.c_double
ps_eisenstein98_zb = np.vectorize(lib.ps_eisenstein98_zb)

args = psargs_t()
args.z       = 0. 
args.h       = 0.7 
args.Omh2    = args_gf.Om0 * args.h**2    
args.Obh2    = 0.04        * args.h**2    
args.ns      = 0.9  
args.sigma8  = 0.8     
args.theta   = 2.725 / 2.7    
args.norm    = 1.    
args.dplus_z = linear_growth(args.z, 0, ctypes.byref(args_gf))       
args.dplus_0 = linear_growth(0.    , 0, ctypes.byref(args_gf))       

lib.init_eisenstein98_zb( ctypes.byref(args) )
# print('z          = ', args.z         )
# print('h          = ', args.h         )
# print('Omh2       = ', args.Omh2      )
# print('Obh2       = ', args.Obh2      )
# print('Onuh2      = ', args.Onuh2     )
# print('Nnu        = ', args.Nnu       )
# print('ns         = ', args.ns        )
# print('sigma8     = ', args.sigma8    )
# print('theta      = ', args.theta     )
# print('dplus_z    = ', args.dplus_z   )
# print('dplus_0    = ', args.dplus_0   )
# print('norm       = ', args.norm      )
# print('include_nu = ', args.include_nu)
# print('z_eq       = ', args.z_eq      )
# print('z_d        = ', args.z_d       )
# print('s          = ', args.s         )
# print('k_silk     = ', args.k_silk    )
# print('param      = ', list(args.param))

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

pktab = np.zeros((1001, 2))
pktab[:, 0] = np.linspace( np.log(1e-04), np.log(1e+04), pktab.shape[0] )
pktab[:, 1] = ps_eisenstein98_zb(pktab[:, 0], ctypes.byref(args))

var8_calc  = variance( np.log(8./args.h), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] )
args.norm  = args.sigma8**2 / var8_calc 
pktab[:,1] = pktab[:,1] + np.log(args.norm)
# print( variance( np.log(8./args.h), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] ), args.sigma8**2 )

# z = np.linspace(0., 10, 11)
# y = linear_growth(z, 0, ctypes.byref(args_gf))
# plt.plot(z, y, 'o-')
# plt.show()

# k = np.logspace(-4, 3, 11)
# y = np.exp( ps_eisenstein98_zb(np.log(k), ctypes.byref(args)) )
# plt.loglog(k, y, '-o')
# plt.show()

# r = np.logspace(-3, 2, 11)
# s = variance( np.log(r), 0, 0, 0, pktab, pktab.shape[0], pktab.shape[1] )
# c = correlation( np.log(r), 0, pktab, pktab.shape[0], pktab.shape[1] )
# plt.loglog(r, s, '-o')
# plt.loglog(r, c, '-o')
# plt.show()
