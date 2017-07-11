#!/usr/bin/env python
"""
This is a wrapper to generate a list of mock catalogs for ImSim
and then run ImSim to generate fits images.
the grid includes surface brightness, stream width (in deg)
and distance of the stream (in kpc). The length of the stream is
set to be 0.2 deg since the LSST CCD size is around 0.3 deg

.txt files are the mock catalog
.info files are the associated information file for each catalog

"""

import os
import numpy as np
import os.path
import doimsim




filedir = 'catalogs'

if not os.path.exists(filedir):
    print 'creating dir', filedir
    os.mkdir(filedir)

sb_vals = [27.0, 29.0] # surface brightness
width_vals = [0.001, 0.01, 0.1] # angular width of the stream in deg
dist_vals = [100., 1000., 10000] # distance in kpc

#set length of the stream be be 0.2 deg (LSST CCD is around 0.3 deg)

for sb in sb_vals:
    for width in width_vals:
        for dist in dist_vals:
            # generate the catalog
            print "running generateInstCat.py to generate the mock catalog ..."
            command = "python generateInstCat.py --surface_brightness=%.2f --angular_width=%f --angular_length=0.2 --distance=%.2f '%s/sb_%i_width_%.4fdeg_d_%ikpc.txt' > %s/sb_%i_width_%.4fdeg_d_%ikpc.info" %(sb, width, dist, filedir, sb, width, dist, filedir, sb, width, dist)
            print command
            os.system(command)
