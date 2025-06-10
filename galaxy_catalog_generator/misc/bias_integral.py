import numpy as np
from scipy.integrate import dblquad
from typing import Protocol, Literal, Any

class _SupportedHaloModel(Protocol):
    def haloMassFunction(self, lnm: float, return_value: Literal["dndlogm", "dndlnm", "dndm", "fs"]) -> float: ...
    def centralCount(self, lnm: float) -> float: ...
    def satelliteCount(self, lnm: float, ) -> float: ...

class _BiasCallable(Protocol):
    def __call__(self, lnm1: float, lnm2: float, *args: Any) -> float: ...

def biasIntegral(
        halomodel: _SupportedHaloModel, 
        biasfunc: _BiasCallable,
        lnma: float,
        lnmb: float = np.log(1e+18),
        bfargs: tuple = (), 
        **kwargs: Any
    ) -> float:
    r"""
    Calculate general galaxy bias integral.

    Parameters
    ----------
    halomodel : HaloModel
        Any object that can be used as a halo model object. It should define these functions of 
        halo mass:

        - ``haloMassFunction(lnm: float, return_value: str) -> float`` for calculating the halo 
            mass-function, were the value ``return_value = 'dndlnm`` would return the value as
            `dn/dlnm`.

        - ``centralCount(lnm: float) -> float`` for calculating central galaxy count

        - ``satelliteCount(lnm: float) -> float`` for calculating satellite galaxy count

    biasfunc : callable 
        Function modelling the inseparable halo-halo bias. In general, it could be a function of 
        the form ``bias(lnm1: float, lnm2: float, *args) -> float``.

        .. note::
            Make sure that the bias model agrees with the halo model settings. 

    lnma : float
        Lower limit of integration.

    lnmb : float, default = log(1e+18)
        Upper linit of integraion. 

    bfargs : tuple, optional
        Additional arguments passed to the bias function.

    **kwargs : any
        Additional keyword arguments passed to the ``scipy.integrate.dblquad`` function.

    Returns
    -------
    bgal2 : float
        Value of the integral - square of the galaxy bias (``b_gal**2``).

    """

    def _integrand(lnm2: float, lnm1: float, *args: Any) -> float:
        
        # calculating halo mass-function...
        dndlnM1 = halomodel.haloMassFunction(lnm1, return_value = "dndlnm")
        dndlnM2 = halomodel.haloMassFunction(lnm2, return_value = "dndlnm")

        # calculating bias values: for the general case, bias is not seperable, but 
        # in some special case, it is separable as `b(m1)*b(m2)`...
        bias = biasfunc(lnm1, lnm2, *args)

        # total halo count...
        galaxyCount1 = halomodel.centralCount(lnm1) + halomodel.satelliteCount(lnm1)
        galaxyCount2 = halomodel.centralCount(lnm2) + halomodel.satelliteCount(lnm2)
        

        retval = dndlnM1 * dndlnM2 * galaxyCount1 * galaxyCount2 * bias
        return retval
    
    # Double integration on the range lnm1, lnm2 = lnma...lnmb:
    kwargs = kwargs | { 
        # default args to the dblquad function:
        "a": lnma, "b": lnmb, "gfun": lnma, "hfun": lnmb, "args": bfargs 
    }
    retval, abserr = dblquad( _integrand, **kwargs )
    return retval # (galaxy bias)^2
