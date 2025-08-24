# Galaxy Catalog Generator

Tools for theoretical calculations related to darkmatter halos.

# ``abacus_halos`` CLI Tool

This CLI tool provides a set of commands for estimating halo mass-function and 
generating galaxy catalogs from halo catalogs of the abacus simulation.


## The ``massfunc`` command

This can be used to estimate halo mass-function from an abacus halo catalog. 
Estimated values will be saved as a CSV table of halo mass ($M_\odot$) and 
mass-function $dn/dm$ ($\rm Mpc^{-3}$). Both values are in natural log format.

### Options:
Option              | Description
--------------------|-----------------------------------------------------------------------
``--siminfo``       | A tuple of a valid abacus simulation name and redshift value (**required**)
``--output-file``   | Filename for the output (**required**)
``--mass-range``    | Values of the left- and right-most mass bins in particle mass unit
``--nbins``         | Number of bins. Must be greater than 2.
``-p``, ``--path``  | Paths to look for simulation catalog files


## The ``halomod`` command

This tool can be used to calculate the optimum values for halo model, based on 
observed values of galaxy density (in $\rm Mpc^{-3}$) and satellite  fraction. 
Other model parameters, including cosmology, are selected corresponding to the 
given abacus simulation details. The output will be saved in YAML format.

### Options:
Option                               | Description
-------------------------------------|-----------------------------------------------------------------------
``--siminfo``                        | A tuple of a valid abacus simulation name and redshift value  (**required**)
``--galaxy-density``, ``--ngal``     | Observed value of the galaxy number density in Mpc^-3  (**required**)
``--satellite-fraction``, ``--fsat`` | Observed value of the satellite galaxy fraction  (**required**)
``--output-file``                    | Filename for the output (**required**)
``--sigma-m``                        | Transition width for central galaxy count function 
``--alpha``                          | Index for the power law satellite count function  
``--scale-shmf``                     | Scale for the SHMF - specify the maximum mass of a subhalo, as a fraction of halo
``--slope-shmf``                     | Slope for the SHMF
``--powerspectrum``, ``--mps``       | Matter power spectrum model used for calculations
``--massfunction``, ``--hmf``        | Halo mass-function model used for calculations or path to the data file
``--mmin-range``                     | Search range for parameter Mmin (or M0) in Msun units
``--m1-range``                       | Search range for parameter M1 in Msun units
``--gridsize``                       | Size of the Mmin-M1 grid used to guess minimum 


## The ``galaxies`` command

This tool can be  used to generate a galaxy catalog using the specified abacus halo 
catalog and halo model parameters.  

**Note**: Always make sure that the halo model parameters match the cosmology setup of 
the simulation used and observation for more realistic results.

### Options:
Option             | Description
-------------------|-----------------------------------------------------------------------
``--siminfo``      | A tuple of a valid abacus simulation name and redshift value (**required**)
``--halomodel``    | Path to the file which store the halo model parameters in YAML format (**required**)
``--output-path``  | Path to the directory to save galaxy catalog files (**required**)
``-p``, ``--path`` | Paths to look for simulation catalog files


### Using the catalog

Galaxy catalogs generated are saved as ASDF files to specified location. These 
files will have two sections 

+ A ``header`` section containing the parameter values used in simulation, and 

+ A ``data`` section containing the catalog data. This table has fields

Attribute           | Type           | Description 
--------------------|----------------|-------------------------------------------
``parentHaloID``    | ``int64``      | ID for the parent halo in the halo catalog
``galaxyPosition``  | ``float64[3]`` | Position coordinates in ${\rm Mpc}/h$
``galaxyMass``      | ``float64``    | Mass in $M_\odot/h$   
``galaxyType``      | ``uint8``      | Galaxy type (`1`: central, `2`: satellite)      