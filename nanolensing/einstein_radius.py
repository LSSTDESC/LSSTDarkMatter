import numpy as np
import scipy.special
import pylab

import cosmology
import cosmology_constants

pylab.ion()

##########

z_L = 0.5
z_S = 2.0

#D_L = cosmology.D_A(z_L) # cm
D_L = 50. * cosmology_constants.PC_TO_CM # cm
D_S = cosmology.D_A(z_S) # cm

##########

def getRhoVir(z):
    # Eq. 6
    return (18. * np.pi**2 + (82 * (getOmegaM(z) - 1.)) - (39. * (getOmegaM(z) - 1.)**2)) * getRhoCrit(z) # M_Earth pc^-3

def getOmegaM(z):
    # Eq. 7
    numerator = cosmology_constants.OMEGA_M * (1. + z)**3
    denominator = cosmology_constants.OMEGA_M * (1. + z)**3 + 1. - cosmology_constants.OMEGA_M
    return numerator / denominator

def getRhoCrit(z):
    # Eq. 8
    c = cosmology_constants.OMEGA_M * (1. + z)**3 + 1. - cosmology_constants.OMEGA_M
    return 0.0924 * cosmology_constants.LITTLE_H**2 * c # M_Earth pc^-3

def getConcentration(m_vir):
    # Eq. 28, M_Vir in solar masses
    return 94. * (m_vir / 1.e6)**(-0.067)

def getFactor(m_vir, gamma, z_vir=0.):
    # Eq. 27.
    c = getConcentration(m_vir)
    a = 3. - gamma
    b = gamma - 2.
    z = c * (gamma - 2.)
    print 'BETA', scipy.special.betainc(a, b, z)
    G = 3.57 * (gamma - 2.) / scipy.special.betainc(a, b, z) # Note that the equation had some weird (-1)**gamma factor, which I set to 1.
    return 810. * (c / 94.) * (m_vir / 1.e6)**(2. / 3.) * (getRhoVir(z_vir) / 4.6)**(1. / 3.) * G

"""
def getSigmaAlpha(m_vir, gamma, D_L, D_S, z_vir=0.):
    # Eq. 22, M_Vir in solar masses
    G = scipy.special.gamma(0.5 * (gamma - 1.)) / (2. * (3. - gamma) * scipy.special.gamma(0.5 * gamma))
    D = 1. - (D_L / D_S)
    F = getFactor(m_vir, gamma, z_vir)
    return 0.88 * G * D * F
"""

def getThetaAlpha(m_01, gamma, D_L, D_S):
    # Eq. 24, m_01 in solar masses, returns angle in microarcseconds
    G = scipy.special.gamma(0.5 * (gamma - 1.)) / (2. * (3. - gamma) * scipy.special.gamma(0.5 * gamma))
    D = 1. - (D_L / D_S)
    #print 'G', G
    #print 'D', D
    return 8.8 * G * D * ((3. - gamma) / (4. * np.pi)) * m_01

def getThetaE(m_01, gamma, D_L, D_S):
    # Eq. 23, returns angle in microarcseconds
    theta_alpha = getThetaAlpha(m_01, gamma, D_L, D_S)
    return (theta_alpha * (0.1 * cosmology_constants.PC_TO_CM / D_L)**(2. - gamma))**(1. / (gamma - 1.))

##########

z_array = np.linspace(0., 5., 1000)

rho_crit_array = getRhoCrit(z_array)
rho_vir_array = getRhoVir(z_array)
omega_m_array = getOmegaM(z_array)

pylab.figure()
pylab.plot(z_array, omega_m_array)

#print D_L / cosmology_constants.MPC_TO_CM
#print D_S / cosmology_constants.MPC_TO_CM

m_vir = 1.e6
m_01 = 1.
gamma = 1.1
#print getConcentration(m_vir)
#print getFactor(m_vir, gamma)
theta_alpha = getThetaAlpha(m_01, gamma, D_L, D_S)
theta_E = getThetaE(m_01, gamma, D_L, D_S)



print theta_alpha
print theta_E
print theta_E**(gamma - 1.) * (0.1 * cosmology_constants.PC_TO_CM / D_L)**(2. - gamma)
