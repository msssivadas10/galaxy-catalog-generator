# Galaxy Catalog Generator

A simple galaxy catalog generator based on a 5 or 3-parameter halo occupation distribution (HOD) model. 

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
``massfunc``.
