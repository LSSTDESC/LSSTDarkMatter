#!/usr/bin/env python
"""
Module to implement likelihood function for quantifying the mass fraction and mass function of dark matter substructure from strongly lensed systems.
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np
import pylab as plt
from scipy import stats
from scipy import integrate
from scipy.integrate import simps,trapz
from scipy.interpolate import interp1d
from scipy.misc import factorial

MMIN=4e6
MMAX=4e9
MLOW=1e8
MHIGH=4e9
P = (MMIN,MMAX,MLOW,MHIGH)

#SIGMA=1e7
SIGMA = 0
ALPHA=1.9
FRAC=0.02
MHALO=1e11
NSTEPS=500
FRAC=1e-2

def mhalo(R):
    return MHALO

def create_mass_array(p=P,nsteps=NSTEPS):
    mmin,mmax,mlow,mhigh = p
    nsteps = nsteps
    m = np.logspace(np.log10(mmin),np.log10(mmax),500)
    mp = np.logspace(np.log10(mlow),np.log10(mhigh),300)
    mm,mmp = np.meshgrid(m,mp)
    return m,mp,mm,mmp

def mu0(alpha, frac, R=1, p=P):
    """ Expectation from the full true mass function (Eqn. 5).
    """
    mmin,mmax,mlow,mhigh=p
    alpha = np.atleast_1d(alpha)
    integral = ( (2-alpha)*(mmax**(1-alpha) - mmin**(1-alpha))) / \
               ( (1-alpha)*(mmax**(2-alpha) - mmin**(2-alpha)))
    integral = np.where(alpha==2,-(mmax**-1 - mmin**-1)/np.log(mmax/mmin),integral)
    integral = np.where(alpha==1,np.log(mmax/mmin)/(mmax - mmin),integral)

    return frac * mhalo(R) * integral

def mu(alpha, frac, R=1, p=P, nsteps=NSTEPS):
    """Expectation value for the number of substructures in an aperature
    with radius R.
    """
    mmin,mmax,mlow,mhigh = p
    m,mp,mm,mmp = create_mass_array(p,nsteps)

    _mu0 = mu0(alpha, frac, R, p)
    if SIGMA == 0:
        _integrand = dP_dm_true(mp,alpha)
        _integral = simps(_integrand,mp)
    else:
        _integrand = dP_dm_conv(mm,mmp,alpha,SIGMA)
        _integral = simps(simps(_integrand,m),mp)
    return _mu0 * _integral

def dP_dm_true(m,alpha=ALPHA,p=P):
    """ Normalized true mass function (Eqn. 6)
    """
    mmin,mmax,mlow,mhigh = p
    m = np.atleast_1d(m)
    ret = ((1-alpha)*m**(-alpha))/(mmax**(1-alpha)-mmin**(1-alpha))
    ret = np.where(alpha==1,(m**-alpha)/np.log(mmax/mmin),ret)
    return np.where(np.isfinite(ret),ret,np.nan)

def dP_dm_conv(m,mp,alpha=ALPHA,sigma=SIGMA):
    """ Convolved mass function (Eqn. 7)
    Parameters:
    ----------
    m  : The range of true masses
    mp : The range of observed masses

    Returns:
    --------
    dP_dm_conv : The convolved probability mass function
    """
    return dP_dm_true(m, alpha)*stats.norm.pdf(m,loc=mp,scale=sigma)

def P_mi(mi, alpha, p=P, nsteps=500):
    """ Probility of measuring mass m_i for a single substructure.
    """
    m,mp,mm,mmp = create_mass_array(p,nsteps)
    if SIGMA == 0:
        top = simps(dP_dm_true(m,alpha),m)
        bottom = simps(dP_dm_true(mp,alpha=alpha),mp)
    else:
        top = simps(dP_dm_conv(m,mi,alpha,sigma=1e7),m)
        bottom = simps(simps(dP_dm_conv(mm,mmp,alpha=alpha,sigma=1e7),m),mp)
    return top/bottom

def prob(m,alpha,p=P,nsteps=NSTEPS):
    if not len(m): return 1
    return np.prod([P_mi(mi,alpha,p,nsteps) for mi in m])

def logprob(m,alpha,p=P,nsteps=NSTEPS):
    if not len(m): return 0
    return np.sum([np.log(P_mi(mi,alpha,p,nsteps)) for mi in m])

def poisson(ns,alpha,frac,R=1,p=(MMIN,MMAX,MLOW,MHIGH),nsteps=NSTEPS):
    """Number of substructures populating a given galactic halo
    fluctuates with a distribution which is Poissonian.

    First factor in the likelihood (Eqn. 1)
    """
    _mu = mu(alpha,frac,R,p,nsteps)
    return stats.poisson.pmf(ns,_mu)

def logpoisson(ns,alpha,frac,R=1,p=(MMIN,MMAX,MLOW,MHIGH),nsteps=NSTEPS):
    _mu = mu(alpha,frac,R,p,nsteps)
    return stats.poisson.logpmf(ns,_mu)

def like(m, alpha, frac, R=1, p=P, nsteps=NSTEPS):
    ns = len(m)
    print ns
    like = poisson(ns,alpha,frac,R,p,nsteps)
    like *= prob(m,alpha,p,nsteps)
    return like

def loglike(m, alpha, frac, R=1, p=P, nsteps=NSTEPS):
    ns = len(m)
    #print ns
    loglike = logpoisson(ns,alpha,frac,R,p,nsteps)
    loglike += logprob(m,alpha,p,nsteps)
    return loglike

def LogPoisson(data, alpha, frac, p=P, nsteps=NSTEPS):
    logprobs = logpoisson(data['nsrc'].flat,alpha,frac)
    return np.sum(logprobs,axis=0)

def LogProb(data, alpha, p=P, nsteps=NSTEPS):
    logprobs =[logprob(m,alpha,p,nsteps) for m in data['src']]
    return np.sum(logprobs,axis=0)

def LogLike(data, alpha, p=P, nsteps=NSTEPS):
    print len(data)
    loglikes =[loglike(m,alpha,p,nsteps) for m in data['src']]
    return np.sum(loglikes,axis=0)

def posterior():
    pass

def sample(num,x,pdf):
    num = int(num)
    cdf = np.cumsum(pdf)
    cdf = np.insert(cdf, 0, 0.)
    cdf /= cdf[-1]

    f = interp1d(cdf, range(0, len(cdf)), bounds_error=False, fill_value=-1)
    index = np.floor(f(np.random.uniform(size=num))).astype(int)
    index = index[index >= 0]
    return x[index]

def simulate(nlens=1,alpha=ALPHA,frac=FRAC, R=1, p=P):
    # First, figure out how many lenses we are sampling
    m,mp,mm,mmp = create_mass_array(p=P,nsteps=NSTEPS)
    pdf = dP_dm_true(m)
    _mu = mu(alpha,frac,R,p)
    lenses = stats.poisson.rvs(_mu,size=nlens)
    out = []
    for i,l in enumerate(lenses):
        masses = sample(l,m,pdf)
        out += [(i,len(masses),masses)]

    names = ['lens','nsrc','src','mass']
    return np.rec.fromrecords(out,names=names)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
