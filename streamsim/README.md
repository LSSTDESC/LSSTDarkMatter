# Description

This is a simple tool to generate mock catalogs for streams and the associated images
  with imSim.

# Requirements:

## Image simulation

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

# To use the code:

* generateInstCat.py

you can use generateInstCat.py to generate mock catalog of streams. use --help to see the input parameters. An example for the terminal command:

 py generateInstCat.py --angular_width 0.05 --angular_length 0.3 --surface_brightness 24 stream0.05_0.3_24.txt

This generates a stream with 0.05 deg x 0.3 deg and surface brightness of 24 mag/arcsec^2, and save it as file stream0.05_0.3_24.txt

The general information of the stream and the simulation will display in the terminal.

Note the code was generating a line of N dwarfs with Gaussian profile to mimic the stream. The half-light radius of the dwarf rh = width/2.35 --0.02 deg in this case and number of dwarf N = length/rh -- 14 in this case.

* doimsim.py

after a mock catalog is generated, you can use doimsim.py to generate images with imSim

* run_generate.py

this code will run generateInstCat.py and doimsim.py to generate a grid of catalogs and the associated images with varius surface brightness, distance, stream width, etc.