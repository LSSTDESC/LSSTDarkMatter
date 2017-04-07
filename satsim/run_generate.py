import os
import numpy as np

abs_mag_vals = [-12.]
r50_vals = [300.]
dist_vals = [5000., 10000., 20000.]

for abs_mag in abs_mag_vals:
    for r50 in r50_vals:
        for dist in dist_vals:
            command = "python generateInstCat.py --absolute_magnitude=%.2f --half_light_radius=%.2f --distance=%.2f 'M_%i_r50_%i_d_%iMpc.txt'" %(abs_mag, r50, dist, abs_mag, r50, dist/1000.)
            print command

            os.system(command)
