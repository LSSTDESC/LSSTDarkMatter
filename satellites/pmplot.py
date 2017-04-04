import down_girardi
import scipy.interpolate
import read_girardi
import numpy as np
import matplotlib.pyplot as plt


def fixer(m_ini, int_imf):
	# remove the points from the interpolated isochrone that are too close in
	# m_ini
	diff = np.diff(int_imf)
	xind = np.where(diff > 1e-10)[0]
	return m_ini[xind], int_imf[xind], xind


def genstar(m_ini, int_imf, dm, N, mags=None, maglim=None, pad=0.5):
	# sample the isochrone
	if maglim is not None:
		if not isinstance(maglim, dict):
			raise Exception('oops')
	else:
		maglim = {}
	m_ini_F, int_imf_F, good_F = fixer(m_ini, int_imf)
	im1, im2 = int_imf_F.min(), int_imf_F.max()
	m1, m2 = m_ini_F.min(), m_ini_F.max()
	II = scipy.interpolate.UnivariateSpline(int_imf_F, m_ini_F, s=0, k=2)

	interps = {}
	# pad = 0.5 # extra padding in magnitude
	mgrid = np.linspace(m1, m2, 100000, True)
	imgrid = np.linspace(im1, im2, 100000)
	xind = np.ones(len(mgrid), dtype=bool)
	for k, v in mags.items():
		magII = scipy.interpolate.UnivariateSpline(m_ini_F, v[good_F], s=0, k=2,
												   ext=2)
		interps[k] = magII
		if k in maglim:
			curxind = (magII(mgrid) + dm) < (maglim[k] + pad)
			xind = xind & curxind
	if xind.sum() == 0:
		return None
	minmass = mgrid[xind].min()
	
	minim = np.nonzero(II(imgrid) >= minmass)[0][0]
	minim = max(0, minim - 1)
	minim = imgrid[minim]  # starting point in the IMF
	res = {}
	cnt = 0
	cnt2 = 0
	while cnt < N:
		#print (cnt)
		masses = II(np.random.uniform(minim, im2, size=int(N)))
		curres = []
		xind = np.ones(N, dtype=bool)
		for k, v in interps.items():
			if k not in res:
				res[k] = []
			curmag = v(masses) + dm
			#1 / 0
			if k in maglim:
				xind = xind & (curmag < maglim[k])
			res[k].append(curmag)

		for k in interps.keys():
			res[k][-1] = res[k][-1][xind]
		cnt += xind.sum()
		cnt2 += N
		#1 / 0
	if len(res)==0:
		return None
	for k in interps.keys():
		res[k] = np.concatenate(res[k])
	return res


def getnorm(iso, dm, maglim):
	# find the luminosity of a system with on average one star above the limit
	nstars = 1e8
	
	dat0 = genstar(iso['M_ini'], iso['int_IMF'], dm, int(nstars),
				   mags={'r': iso['r'], 'g': iso['g']},	
				   #maglim={'r': maglim}
				   )
	Ms = dat0['r'] - dm
	totlum = -2.5 * np.log10((10 ** (-Ms / 2.5)).sum()/nstars)
	frac = (dat0['r']<maglim).sum() / nstars
	totlum = totlum + 2.5*np.log10(frac)
	return totlum

maglim = 25
cache = {}
iso = down_girardi.getit(10, -2, "SDSS ugriz")


def simSat(iso, mv, dm, maglim):
	# simulate a satellite with luminosity mv at a distance modulus of dm
	# with the magnitude limit of maglim
	if (dm, maglim) not in cache:
		cache[dm, maglim] = getnorm(iso, dm, maglim)
	norm = cache[dm, maglim]

	# expected number of stars above the limit
	Nstars = 10**((norm - mv) / 2.5)
	Nstars1 = np.random.poisson(Nstars)
	# actual number of stars1
	print (Nstars,Nstars1)
	dat = genstar(iso['M_ini'], iso['int_IMF'], dm, Nstars1,
				  mags={'r': iso['r'],
						'g': iso['g']},
				  maglim={'r': maglim})
				  
	return dat


def getVelPrec(mv, dist, maglim, muPmI):
	# args luminosity[mags], distance[kpc]
	dm = 5 * np.log10(dist * 1e3) - 5
	dat = simSat(iso, mv, dm, maglim)
	if dat is None:
		return np.nan
	muerr = 10**muPmI(dat['r'])
	muerr = np.sqrt(1. / (1. / muerr**2).sum())
	rverr = muerr * dist * 4.74
	return rverr

def getVelPrec1(mv, dist, maglim, muPmI):
	# args luminosity[mags], distance[kpc]
	dm = 5 * np.log10(dist * 1e3) - 5
	dat = simSat(iso, mv, dm, maglim)
	if dat is None:
		return np.nan

	muerr = 10**muPmI(dat['r'])
	rverr = muerr * dist * 4.74
	rvsig = 5
	ret = 1/np.sqrt(2) /rvsig / np.sqrt(np.sum(1/(rvsig**2+rverr**2)**2))
	return ret


# LSST precision in proper motions
mag_lsst = np.arange(16, 25.1, 1.)
mu_lsst = np.array([0.15, 0.15, 0.15, 0.15, 0.15,
					0.18, 0.20, 0.4, 1., 3.])  # mas/yr
# pm interpolator
muPmI = scipy.interpolate.UnivariateSpline(
	mag_lsst, np.log10(mu_lsst), s=0, k=2)


def plotter(plotname):
	plt.clf()
	maglim = 25
	lums = [-10, -7.5, -5, -2.5, 0]
	colors = ['blue', 'green', 'red']
	for ii, d in enumerate([30, 100, 300]):
		res = []
		for l in lums:
			print('doing distance %f, Mr %f' % (d, l))
			res.append(getVelPrec1(l, d, maglim, muPmI))
		plt.plot(lums, res, color=colors[ii], label='%d kpc' % d)
	plt.plot([-10,0,],[3,3],':',color='black')
	plt.gca().set_yscale('log')
	plt.xlabel(r'M$_r$ [mag]')
	plt.ylabel(r'$\sigma_V$ [km/s]')
	plt.legend()

	plt.savefig(plotname)
