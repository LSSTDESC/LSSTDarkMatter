from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from ugali.analysis.isochrone import factory as isochrone_factory
from ugali.utils.projector import dist2mod, mod2dist

mag_lim = 26.5
rad_to_arcsec = 206265.

age = 12.
z = 0.0001  # feh = -2.1

dist = np.array([50, 100, 500, 1000, 5000, 10000])  # kpc
m_star = 10**np.arange(3, 7.5, 0.5)  # Msun
r50 = np.logspace(0, 3, 50)  # pc

M, R = np.meshgrid(m_star, r50)

for d in dist:
    print 'd = %.2f' % d

    m_plt = []
    r_plt = []

    # get isochrone
    mu = dist2mod(d)
    iso = isochrone_factory('Padova', age=age, distance_modulus=mu, z=z)

    for m in m_star:
        print 'm = %.2e' % m
        # calculate number of stars above mag lim
        g, r = iso.simulate(m)
        g = g[g < mag_lim]
        r = r[r < mag_lim]
        nstars = len(g) / 2  # use g for now

        # calculate stellar density
        area = 2 * np.pi * (r50/(d*1000) * rad_to_arcsec)**2
        angular_density = nstars / area  # stars / pc^2

        r_det = r50[(angular_density < 1.0) & (nstars > 20)]
        r_plt.extend(r_det)
        m_plt.extend([m] * len(r_det))

        print 'ndwarfs = %i' % len(r_det)

    plt.figure()
    plt.scatter(M, R, facecolors='none', edgecolor='darkslateblue', s=50)
    plt.scatter(m_plt, r_plt, color='darkslateblue', edgecolor='none', s=50)
    plt.xlabel('L (Lsun)')
    plt.ylabel('r50 (pc)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e2, 1e8)
    plt.ylim(r50.min() - 100, r50.max() + 100)
    plt.title('Distance = %.2d kpc' % d)
    plt.savefig('resolvable_d_%.2f.png' %d)
    