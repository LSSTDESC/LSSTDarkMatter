
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

after ugali and dsphsim is installed, you can use generateInstCat.py to generate mock catalog of streams. use --help to see the input parameters. An example for the terminal command:

 py generateInstCat.py --angular_width 0.05 --angular_length 0.3 --surface_brightness 24 stream0.05_0.3_24.txt

This generates a streamw with 0.05 deg x 0.3 deg and surface brightness of 24 mag/arcsec^2, and save it as file stream0.05_0.3_24.txt

You will also see the general information below in the terminal

distance_modulus 17.5
angular_width: 0.05 deg
angular_length: 0.3 deg
angular_radius: 0.0212765957447 deg
going to simulate 14 dwarfs
stream extension: 194400.0 in arcsec^2
surface brightness 24.0
total apparent magnitude 10.7782593485
total absolute magnitude -6.72174065148
single absolute magnitude -3.85642056228
center ra, dec: 54.0 -54.0
each dwarf has 13442 stars
the stream has 188188 stars

Note the code was generating a line of N dwarfs with Gaussian profile to mimic the stream. The half-light radius of the dwarf rh = width/2.35 --0.02 deg in this case and number of dwarf N = length/rh -- 14 in this case.
