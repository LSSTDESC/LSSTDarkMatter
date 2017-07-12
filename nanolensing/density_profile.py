import numpy as np
import pylab

import cosmology_constants

pylab.ion()

##########

"""
x = np.linspace(1.e-3, 3., 1000)
rho = 1. / (x * (1 + x)**2) 

pylab.figure()
pylab.xscale('log')
pylab.yscale('log')
pylab.plot(x, rho)
"""

def profileParameters(m_200):
    """
    m in solar masses
    should generalize the slope
    """
    rho_0 = np.sqrt(m_200 * cosmology_constants.M_SOLAR_TO_G * (3. / (8. * np.pi)) * 200. * cosmology_constants.RHO_CRIT) # g cm^(-3/2) 
    r_max = (rho_0 / (200. * cosmology_constants.RHO_CRIT))**(2. / 3.)
    
    r_max *= cosmology_constants.PC_TO_CM**(-1.) # pc
    rho_0 *= cosmology_constants.M_SOLAR_TO_G**(-1.) * cosmology_constants.PC_TO_CM**(3. / 2.) # M_sol 

    return rho_0, r_max

##########

m_200 = 1.
rho_0, r_max = profileParameters(m_200)

#print (8. * np.pi / 3.) * rho_0 * r_max**(3. / 2.) * cosmology_constants.M_SOLAR_TO_G**(-1.)

#print r_max * cosmology_constants.PC_TO_CM**(-1.)

r = np.linspace(1.e-8, 1., 100000) # pc
rho =  rho_0 * r**(-1.5) # M_sol pc^-3

pylab.figure()
pylab.xscale('log')
pylab.yscale('log')
pylab.plot(r, rho)
pylab.axhline(200 * cosmology_constants.RHO_CRIT * cosmology_constants.M_SOLAR_TO_G**(-1) * cosmology_constants.PC_TO_CM**3, c='red')

"""

# A different strategy

bins = np.linspace(-20, 20, 400)
xx, yy, zz = np.meshgrid(bins, bins, bins)
rr = np.sqrt(xx**2 + yy**2 + zz**2)
rho = rr**(-1.5) * (rr < 5.)
sigma = np.sum(rho, axis=2)

x, y = np.meshgrid(bins, bins)
#r = np.sqrt(xx[0,:,:]**2 + yy[:,0,:]**2)
r = np.sqrt(x**2 + y**2)

pylab.figure()
#pylab.pcolor(xx[0,:,:], yy[:,0,:], sigma, cmap='binary')
pylab.pcolor(x, y, sigma, cmap='binary')



pylab.figure()
pylab.xscale('log')
pylab.yscale('log')
pylab.scatter(r.flatten(), sigma.flatten(), c='blue', edgecolor='none')

r_3d = np.linspace(1.e-2, 1.e2, 1000)
density = 70. * r_3d**(-0.5)

pylab.plot(r_3d, density, c='red')
"""
