
import os, logging, numpy as np
from dataclasses import dataclass
from typing import Literal
from numpy.typing import NDArray
from Corrfunc.theory.DD import DD

@dataclass(frozen = True, slots = True)
class PairCountData:
    r"""
    Object storing pair count data.
    """

    ND1  : int
    ND2  : int
    NR1  : int
    NR2  : int
    bins : NDArray[np.float64]
    D1D2 : NDArray[np.uint64]
    D1R2 : NDArray[np.uint64]
    D2R1 : NDArray[np.uint64]
    R1R2 : NDArray[np.uint64]

    @classmethod
    def countPairs(
            cls,
            D1       : NDArray[np.float64],
            D2       : NDArray[np.float64], 
            rbins    : list[float], 
            boxsize  : float,
            periodic : bool  = False, 
            diff_r   : bool  = False,
            rsizefrac: float = 3., 
            nthreads : int   = 0,
            seed     : int   = None,
        ):
        r"""
        Count the number of pairs from catalogs `D1` and `D2`.

        Parameters
        ----------
        D1, D2 : ndarray of shape (N, 3)
            Catalogs. If `D2` is None, auto correlation for `D1` catalog is calculated. Catalog positions
            should be in the range `[-boxsize/2, boxsize/2]`.

        rbins : sequence of float
            Distance bin edges.
        
        boxsize : float
            Size of the bounding box (cube) containing all points. 
        
        periodic : bool, default=False
            Tell if to use periodic wrapping of space.

        diff_r : bool, default=False
            Tell if to use two different random catalogs corresponding to the catalogs `D1` and `D2`. It 
            should be true, if the data catalogs are from different samples. 

        rsizefrac : int, default=3
            Size of the random catalog is calculated as this factor times the data size.

        nthreads : int, optional
            Number of threads to use. If set to 0 or not specified, use all available threads.

        seed : int, optional
            Random number seed value.

        Returns
        -------
        obj : `PairCountData`
            Object containing the pair counts. This can be used to calcalute correlation function.

        Usage
        ----- 
        >>> import numpy as np
        >>> # generating random points (correlation = 0)
        >>> boxsize = 100
        >>> D1 = np.random.uniform(-boxsize/2, boxsize/2, (10000, 3))
        >>> D2 = np.random.uniform(-boxsize/2, boxsize/2, (10000, 3))
        >>> rbins = [5, 10, 15]
        >>> pc = PairCountData.countPairs(D1, D2, rbins, boxsize, False)

        Correlation can be calculated as (using LS estimator)

        >>> xi, err = pc.correlation()
        >>> xi
        array([-0.00106298,  0.00071791])

        This is close to the expected value for random points, 0.

        """

        logger = logging.getLogger()

        D1 = np.asarray(D1) 
        D2 = np.asarray(D2) if D2 is not None else D1
        assert D1.ndim == 2 and D1.shape[1] > 2, "D1 array should be a 2d array of 3 columns"
        assert D2.ndim == 2 and D2.shape[1] > 2, "D2 array should be a 2d array of 3 columns"

        rng = np.random.default_rng(seed)
        if diff_r:
            # Using two different different random catalogs
            
            Nr1 = int( rsizefrac*D1.shape[0] )
            logger.info(f"generating random catalog R1 of size {Nr1}...")
            R1  = rng.uniform(-0.5*boxsize, 0.5*boxsize, size = [Nr1, 3])

            Nr2 = int( rsizefrac*D2.shape[0] )
            logger.info(f"generating random catalog R2 of size {Nr2}...")
            R2  = rng.uniform(-0.5*boxsize, 0.5*boxsize, size = [Nr2, 3])
        else:
            # Using same random catalog

            Nr1 = int( rsizefrac*max( D1.shape[0], D2.shape[0] ) )
            logger.info(f"generating random catalog R1 of size {Nr1}...")
            R1  = rng.uniform(-0.5*boxsize, 0.5*boxsize, size = [Nr1, 3])
            R2  = R1

        # Common args to all pair count calls:
        kwargs = dict(
            nthreads = nthreads or os.cpu_count(), 
            binfile  = rbins, 
            periodic = periodic, 
            boxsize  = boxsize,
            verbose  = False,
        ) 

        # -- Calculating pair counts D1D2, D1R2, D2R1 and R1R2:

        # Calculating number of pairs between D1 and R2 catalogs, D1R2, as it is common if D2 is
        # same as D1 or not...  
        logger.info(f"counting pairs of D1 and R2...")
        D1R2 = DD(
            autocorr = 0, 
            X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
            X2 = R1[:,0], Y2 = R1[:,1], Z2 = R1[:,2], 
            **kwargs, 
        )
        D1R2 = D1R2["npairs"]

        if D2 is D1:
            # Here, D2 catalog is same as D1 catalog, so D1D2 is the number of pairs in the D1
            # catalog only (``autocorr = 1``) and D2R1 is same as D1R2... 
            logger.info(f"counting pairs of D1 and D1...")
            D1D2 = DD(
                autocorr = 1, 
                X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
                **kwargs, 
            )
            D1D2 = D1D2["npairs"]

            D2R1 = D1R2
        else:
            # D1 and D2 are different, calculating D1D2 and D2R1... 
            logger.info(f"counting pairs of D1 and D2...")
            D1D2 = DD(
                autocorr = 0, 
                X1 = D1[:,0], Y1 = D1[:,1], Z1 = D1[:,2], 
                X2 = D2[:,0], Y2 = D2[:,1], Z2 = D2[:,2], 
                **kwargs, 
            )     
            D1D2 = D1D2["npairs"]
        
            logger.info(f"counting pairs of D2 and R1...")
            D2R1 = DD(
                autocorr = 0, 
                X1 = R1[:,0], Y1 = R1[:,1], Z1 = R1[:,2], 
                X2 = D2[:,0], Y2 = D2[:,1], Z2 = D2[:,2], 
                **kwargs, 
            )
            D2R1 = D2R1["npairs"]

        logger.info(f"counting pairs of R1 and R2...")
        if diff_r:
            R1R2 = DD(
                autocorr = 0, 
                X1 = R1[:,0], Y1 = R1[:,1], Z1 = R1[:,2], 
                X2 = R2[:,0], Y2 = R2[:,1], Z2 = R2[:,2], 
                **kwargs, 
            )
            R1R2 = R1R2["npairs"]
        else:
            R1R2 = DD(
                autocorr = 1, 
                X1 = R1[:,0], Y1 = R1[:,1], Z1 = R1[:,2], 
                **kwargs, 
            )
            R1R2 = R1R2["npairs"]

        obj = cls(
            ND1  = D1.shape[0], 
            ND2  = D2.shape[0], 
            NR1  = R1.shape[0], 
            NR2  = R2.shape[0],
            bins = np.asarray(rbins),
            D1D2 = D1D2, 
            D1R2 = D1R2, 
            D2R1 = D2R1, 
            R1R2 = R1R2,
        )
        return obj    

    def correlationNatural(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        r"""
        Calculating 2-point correlation function based on Natural estimator.

        .. math::

            \xi = \frac{D_1D_2}{R_1R_2}  - 1 

        Returns
        -------
        xcf : array_like
            Correlation function values.

        err : array_like
            Error in the estimated values, assuming Poisson noise in counts.
        
        """

        logger = logging.getLogger()
        logger.info(f"calculating correlation function using 'natural' estimator...")

        # normalizing counts:
        d1d2 = self.D1D2 / ( self.ND1*self.ND2 )
        r1r2 = self.R1R2 / ( self.NR1*self.NR2 )
        
        # correlation function:
        xcf  = d1d2 / r1r2 - 1.

        # relative error, assuming poisson noise:
        err  = 1 / np.sqrt(self.D1D2)  + 1 / np.sqrt(self.R1R2)
        
        return xcf, np.abs(xcf + 1.) * err

    def correlationDavisPeebles(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        r"""
        Calculating 2-point correlation function based on Davis-Peebles estimator.

        .. math::

            \xi = \frac{D_1D_2}{\frac{1}{2}(D_1R_2 + D_2R_1)}  - 1 

        Returns
        -------
        xcf : array_like
            Correlation function values.

        err : array_like
            Error in the estimated values, assuming Poisson noise in counts.
        
        """
        
        logger = logging.getLogger()
        logger.info(f"calculating correlation function using 'Davis-Peebles' estimator...")

        # normalizing counts:
        d1d2 = self.D1D2 / ( self.ND1*self.ND2 )
        d1r2 = self.D1R2 / ( self.ND1*self.NR2 )
        d2r1 = self.D2R1 / ( self.ND2*self.NR1 )
        
        # correlation function:
        xcf  = 2*d1d2 / (d1r2 + d2r1) - 1.
        
        # relative error, assuming poisson noise:
        err  = (
            1 / np.sqrt(self.D1D2) 
                + ( np.sqrt(self.D1R2) + np.sqrt(self.D2R1) ) / ( self.D1R2 + self.D2R1 )
        )  
        
        return xcf, np.abs(xcf + 1.) * err

    def correlationHamilton(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        r"""
        Calculating 2-point correlation function based on Hamilton estimator.

        .. math::

            \xi = \frac{D_1D_2}{D_1R_2} \frac{R_1R_2}{D_2R_1} - 1 

        Returns
        -------
        xcf : array_like
            Correlation function values.

        err : array_like
            Error in the estimated values, assuming Poisson noise in counts.
        
        """
        
        logger = logging.getLogger()
        logger.info(f"calculating correlation function using 'Hamilton' estimator...")

        # normalizing counts:
        d1d2 = self.D1D2 / ( self.ND1*self.ND2 )
        d1r2 = self.D1R2 / ( self.ND1*self.NR2 )
        d2r1 = self.D2R1 / ( self.ND2*self.NR1 )
        r1r2 = self.R1R2 / ( self.NR1*self.NR2 )
        
        # correlation function:
        xcf  = (d1d2 / d1r2) * (r1r2 / d2r1) - 1.
        
        # relative error, assuming poisson noise:
        err  = ( 
            1 / np.sqrt(self.D1D2) + 1 / np.sqrt(self.D1R2) 
                                   + 1 / np.sqrt(self.R1R2) 
                                   + 1 / np.sqrt(self.D2R1)
        ) 
        
        return xcf, np.abs(xcf + 1.) * err

    def correlationLandySzalay(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        r"""
        Calculating 2-point correlation function based on Landy-Szalay estimator.

        .. math::

            \xi = \frac{D_1D_2 - D_1R_2 - D_2R_1 + R_1R_2}{R_1R_2}

        Returns
        -------
        xcf : array_like
            Correlation function values.

        err : array_like
            Error in the estimated values, assuming Poisson noise in counts.
        
        """
        
        logger = logging.getLogger()
        logger.info(f"calculating correlation function using 'Landy-Szalay' estimator...")

        # normalizing counts:
        d1d2 = self.D1D2 / ( self.ND1*self.ND2 )
        d1r2 = self.D1R2 / ( self.ND1*self.NR2 )
        d2r1 = self.D2R1 / ( self.ND2*self.NR1 )
        r1r2 = self.R1R2 / ( self.NR1*self.NR2 )
        
        # correlation function:
        y    = (d1d2 - d1r2 - d2r1)
        xcf  = y / r1r2 + 1
        
        # relative error, assuming poisson noise:
        err  = (( 
                (1 / ( self.ND1*self.ND2 )) / np.sqrt(self.D1D2) 
                    + (1 / ( self.ND1*self.NR2 )) / np.sqrt(self.D1R2) 
                    + (1 / ( self.ND2*self.NR1 )) / np.sqrt(self.D2R1) 
            ) / y 
            + 1 / np.sqrt(self.R1R2)
        )
        
        return xcf, np.abs(xcf) * err

    def correlation(
            self, 
            estimator: Literal["ls", "nat", "dp", "ham"] = "ls",
        ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        r"""
        Calculating 2-point correlation function based on the specified estimator.

        Parameters
        ----------
        estimator : { 'ls', 'nat', 'dp', 'ham' }, optional
            Specify the estimator. Available estimators are Landy-Szalay (`ls`, default), natural (`nat`), 
            Davis-Peebles (`dp`) and Hamilton (`ham`).

        Returns
        -------
        xcf : array_like
            Correlation function values.

        err : array_like
            Error in the estimated values, assuming Poisson noise in counts.
        
        """
        
        estimator = estimator.lower()
        if estimator in [ "nat", "natural" ]:
            return self.correlationNatural()
        if estimator in [ "dp", "davis-peebles" ]:
            return self.correlationDavisPeebles()
        if estimator in [ "ham", "hamilton" ]:
            return self.correlationHamilton()
        if estimator in [ "ls", "landy-szalay" ]:
            return self.correlationLandySzalay()
        
        raise ValueError(f"unknown correlation function estimator: {estimator}")
        