# Galaxy Catalog Generator

A simple galaxy catalog generator based on a 5 or 3-parameter halo occupation 
distribution (HOD) model. This model uses a smooth step function for modelling 
the average number of central galaxies in a halo of mass $M$, 

$$
    N_{\rm cen}(M) = 0.5 \left[ 
            1 + {\rm erf}\left( 
                \frac{ \ln M - \ln M_{\rm min}}{ \sigma_M } 
        \right) 
    \right]
$$

For $\sigma_M = 0$, this becomes the usual step function. Average satellite count 
is given as a power law 

$$
    N_{\rm sat} = \left( \frac{ M - M_0}{ M_1 } \right)^\alpha 
$$

Mass distribution of the satellite galaxies is given by the subhalo mass-function, 
which is taken as a bounded Pareto distribution

$$
    n(m | M) \propto \left( \frac{m}{\beta M} \right)^{-\gamma} 
$$

where the scale parameter $\beta$ specify the maximum mass of the satellite galaxy 
and $\gamma$ is the logarithmic slope of the function. Usually $\beta=0.5$ and $\gamma=2$ are used.   

## Using `abacus_galaxies.py`

To generate galaxy catalogs from abacus halo catalogs, create a parameter file 
containing the halo model parameters and other settings in a valid YAML format. 
Then run the script as

```sh
$ python3 abacus_galaxies.py parameter_file.yml
```

Required fields in the parameter file are 

+ ``Mmin``: Specify the minimum halo mass (in Msun) required to host a galaxy
+ ``sigmaM``: Transition width for central galaxy count relation. 0 value 
    is results in a step function.  
+ ``M0``: Minimum mass (in Msun) for having non-zero satellites. It can be taken
    the same as ``Mmin``.
+ ``M1``: Amplitude of the satellite count function.
+ ``alpha``: Slope of the satellite count function. Usually 1. 
+ ``scaleSHMF``: Specify the maximum satellite mass as fraction of halo mass. 
+ ``slopeSHMF``: Slope of the sub-halo mass-function power law.
+ ``catpath``: Path to the halo catalogs. 
+ ``outpath``: Path to the output location to save.
+ ``powerspec``: Power spectrum model to use (default is ``eisenstein98_zb``)
+ ``massfunc``: Halo mass-function model to use (default is ``tinker08``)

To match with observation, ``Mmin`` and ``M1`` are calculated so that the galaxy 
density and satellite fraction values calculated using these values matches with 
the observed values. To do this ``hod_optimizer.py`` script can be used as

```sh
$ python3 hod_optmizer.py --OPTION1=VALUE1 --OPTION2=VALUE2 ...
```

Important options are:

+ ``simname``: Name of the abacus simulation (to get cosmology parameters).
+ ``redshift``: Redshift of the catalog.
+ ``galdens``: Observed value of galaxy density in ${\rm Mpc}^{-3}$.
+ ``satfrac``: Observed value of satellite fraction.

Other parameters are ``sigma`` (same as ``sigmaM``), ``alpha``, ``powerspec`` and 
``massfunc``. Running this will print the calculated optimum parameters to the 
stdout in a YAML compatible format. This can be copied to the parameter file for 
using with the catalog generator. Also, the values are saved as the ASDF file 
``hod_optimizer_result.asdf``.

## Using the catalog

Galaxy catalogs generated are saved as ASDF files to specified location. These 
files will have a ``header`` section containing the parameter values and a 
``data`` section containing the catalog data. This table has fields

+ ``parentHaloID``: an ``int64`` ID for the parent halo in the halo catalog. 
+ ``galaxyPosition``: galaxy position coordinates in ${\rm Mpc}/h$ units, as 
    ``float64[3]`` vector.
+ ``galaxyMass``: galaxy mass in $M_\odot/h$ units, as ``float64``. 
+ ``galaxyType``: a ``uint8`` value indicating the type of the galaxy. Central 
    galaxy have value `1` and satellites have `2`.

Another ``statistics`` field in the file store the number of halos uesd, central 
and satellite galaxies in the file, which may be useful later. 