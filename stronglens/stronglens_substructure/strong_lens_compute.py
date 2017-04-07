import numpy as np
import matplotlib.pyplot as plt
import re_me

da=np.genfromtxt('slacs_lenses.txt')
ind=np.arange(0, len(da[0, :]), 2)
Re=da[0, ind]
zl=da[1, ind]
zs=da[1, ind+1]
sigma_a=da[2, ind]

ff = open('strong_lens_compute.txt','w')
ff.write('## Computing properties for the strong lensing systems in Bolton 2006, Apj 638:703-724\n')
ff.write('## ZS, ZL, sigma_a (km/s), theta_E (asec), sigma_c (solar mass/Mpc^2), m_einstein (solar mass), ring mass for dr=0.6" (solar mass)\n')
for ii in range(len(zs)):
    
    DS, DL, DLS=re_me.calc_ds(zs[ii], zl[ii])
    theta_E=re_me.theta_sigma(DLS, DS, sigma_a[ii])
    mE=re_me.Me_theta(theta_E, DL, DS, DLS)
    sigma_c=re_me.Sigma_cden(DL, DS, DLS)

    dtheta=0.6
    r_E=re_me.r_fromtheta(theta_E, DL)
    dr=re_me.r_fromtheta(dtheta, DL)
    ring_mass=re_me.ring_mass(r_E, dr, sigma_c)
    ff.write( '%f %f %f %f %e %e %e \n'%(zs[ii], zl[ii], sigma_a[ii], theta_E, sigma_c, mE, ring_mass) )

ff.close()
