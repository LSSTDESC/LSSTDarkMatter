# Simulating Dwarf Galaxies

The purpose of this module is to explore image simulation and analysis of dwarf galaxies. It includes simulation of both the resolved stars and unresolved light associated with the galaxies.

# Usage

The easiest way to run this code is on NERSC where all the software is pre-installed. The work flow looks something like this:
```
> ssh <username>@cori.nersc.gov
> git clone https://github.com/LSSTDESC/LSSTDarkMatter.git
> cd LSSTDarkMatter/satsim
> source setup-nersc.sh
> python generateInstCat.py instcat.txt
> python doimsim.py instcat_g.txt
> python doimsim.py instcat_r.txt
```
This should create a directory called `fits` that contains a set of CCD images. The dwarf galaxy is centered in CCD sensor (1,1) (`lsst_e_1_R_2_2_S_1_1_g.fits`).

# Requirements:

## Image simulation

* Installation instructions
https://confluence.lsstcorp.org/display/SIM/Catalogs+and+MAF

* Main LSST installation 
eups distrib install -t v13_0 lsst_distrib

* LSST sims installation
eups distrib install -t sims_2_3_4  lsst_sims

* imsim
https://github.com/LSSTDESC/imSim

## Catalog generation

* ugali
https://github.com/DarkEnergySurvey/ugali

* dsphsim
https://github.com/kadrlica/dsphsim