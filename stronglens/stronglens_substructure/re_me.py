import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatwCDM
import astropy.cosmology
from astropy import units as u
cosmo = FlatwCDM(H0=70, Om0=0.3)

def theta_e(M, DL, DS, DLS):
    G = 4.519e-48 #Mpc^3/s^2/M_sun
    m200 = M*1e14
    cc=9.71561*10.0**(-15)# Mpc/sec

    theta_2=4.0*G*m200/cc**2*DLS/DL/DS
    theta=np.sqrt(theta_2)*3600.0*180.0/np.pi
    return theta

def Me_theta(theta, DL, DS, DLS):
    G = 4.519e-48 #Mpc^3/s^2/M_sun
    cc=9.71561*10.0**(-15)# Mpc/sec

    theta_2=(theta*np.pi/3600.0/180.0)**2
    me=theta_2/4.0/G*cc**2/DLS*DL*DS

    return me

def theta_sigma(DLS, DS, sigma_a):
    c=3.0*10.0**5
    theta=4.0*np.pi*(sigma_a/c)**2*(DLS/DS)
    theta=theta*3600.0*180.0/np.pi
    return theta

def calc_ds(zs, zl):
    DS=cosmo.comoving_distance(zs).value/(1.0+zs)
    DL=cosmo.comoving_distance(zl).value/(1.0+zl)
    DLS=1.0/(1.0+zs)*( DS*(1.0+zs)-DL*(1.0+zl) )
    return DS, DL, DLS

def Sigma_cden(DL, DS, DLS):
    G = 4.519e-48 #Mpc^3/s^2/M_sun
    cc=9.71561*10.0**(-15)# Mpc/sec
		    
    sigma_cden=cc**2/4.0/np.pi/G*(DS/DL/DLS)
    return sigma_cden

def ring_mass(theta, dtheta, sigma_cden):
    sigma_dm=0.5*0.63*sigma_cden
    mdm=4.0*np.pi*theta*sigma_dm*dtheta
    return mdm

def r_fromtheta(theta, D):
    r=theta/3600.0/180.0*np.pi*D
    return r

def theta_fromr(r, D):
    theta=r/D/np.pi*180.0*3600.0
    return theta
