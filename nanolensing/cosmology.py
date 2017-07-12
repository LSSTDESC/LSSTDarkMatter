import numpy as np
import scipy.interpolate

import cosmology_constants

############################################################

def E(z):
    return np.sqrt(cosmology_constants.OMEGA_M * (1. + z)**3 + cosmology_constants.OMEGA_LAMBDA)

############################################################

D_H = cosmology_constants.C_LIGHT / cosmology_constants.HUBBLE # cm

############################################################

z_array = np.linspace(0., 100., 1000000)
dz = z_array[1] - z_array[0]
D_C = scipy.interpolate.interp1d(z_array, 
                                 D_H * dz * np.cumsum(E(z_array)**(-1)))

############################################################

def D_L(z):
    return (1. + z) * D_C(z)

############################################################

def D_A(z_1, z_2=None):
    """
    https://arxiv.org/pdf/astro-ph/9905116.pdf
    Equations 18 and 19
    """
    if z_2 is not None:
        return (1. + z_2)**(-1) * (D_C(z_2) - D_C(z_1))
    else:
        return (1. + z_1)**(-1) * D_C(z_1)        

############################################################

if __name__ == '__main__':
    
    import pylab

    pylab.ion()

    z_array = np.linspace(0., 10., 10000)

    pylab.figure()
    pylab.plot(z_array, D_C(z_array) / D_H)
    pylab.plot(z_array, D_A(z_array) / D_H)
    pylab.plot(z_array, D_L(z_array) / D_H)

############################################################
