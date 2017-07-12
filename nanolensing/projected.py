import numpy as np
import scipy.integrate
import pylab

pylab.ion()

#####

def integrate(big_r, alpha):
    l = np.linspace(-50 * big_r, 50 * big_r, 10000)
    y = (l**2 + big_r**2)**(alpha / 2.)
    return scipy.integrate.simps(y, l)

#####

alpha = -1.5
big_r_array = np.arange(0.1, 1, 0.1)
sigma_array = np.empty(len(big_r_array))
for ii in range(0, len(big_r_array)):
    sigma_array[ii] = integrate(big_r_array[ii], alpha)
    

pylab.figure()
pylab.xscale('log')
pylab.yscale('log')
pylab.plot(big_r_array, sigma_array)

print (np.log(sigma_array[1]) - np.log(sigma_array[0])) / (np.log(big_r_array[1]) - np.log(big_r_array[0])) 
