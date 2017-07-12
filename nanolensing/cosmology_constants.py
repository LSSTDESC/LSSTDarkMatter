"""
Useful constants
"""

import numpy as np

# Unit conversions
MPC_TO_CM = 3.086e24
KM_TO_CM = 1.e5 
YR_TO_S = 3.154e7
SR_TO_DEG2 = (180. / np.pi)**2
M_SOLAR_TO_G = 1.988e33 
PC_TO_CM = 3.086e18

# Physical constants
C_LIGHT = 2.998e10 # cm s^-1
G_NEWTON = 6.674e-8 # cm^3 g^-1 s^-2

# Cosmological parameters
LITTLE_H = 0.7
HUBBLE = LITTLE_H * 100. * KM_TO_CM / MPC_TO_CM # (km s^-1 Mpc^-1) to (s^-1) 
OMEGA_M = 0.3
OMEGA_LAMBDA = 0.7
RHO_CRIT = 3. * HUBBLE**2 / (8. * np.pi * G_NEWTON) # g cm^-3
