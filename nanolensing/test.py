import numpy as np
import scipy.interpolate
import pylab

pylab.ion()

##########

#s = 128 # 2**7
s = 2**8

n, m = np.meshgrid(np.fft.fftshift(np.fft.fftfreq(s)) * s, np.fft.fftshift(np.fft.fftfreq(s)) * s)

kappa = np.zeros([s, s])
kappa[int(s / 2.), int(s / 2.)] = 100.
kappa[int(s / 4.), int(s / 4.)] = 100.

#r = np.sqrt((n - 0.5)**2 + (m - 0.5)**2)
#kappa = 1. * r**(-0.5)

pylab.figure()
#pylab.pcolor(2. * (np.cos(m * np.pi / s) + np.cos(n * np.pi / s) - 2.), cmap='binary')
#pylab.pcolor(kappa, cmap='binary')
#pylab.contourf(n, m, kappa)
pylab.contourf(n, m, kappa, 100, cmap='binary')
pylab.colorbar()
#pylab.xlim(0., s)
#pylab.ylim(0., s)
pylab.xlim(-s / 2, s / 2)
pylab.ylim(-s / 2, s / 2)

kappa_fft = np.fft.fftshift(np.fft.fft2(kappa))
k_squared = 2. * (np.sin(m * np.pi / s)**2 + np.sin(n * np.pi / s)**2)
k_squared[k_squared == 0.] = 1.
phi_fft = kappa_fft / k_squared
#phi_fft = kappa_fft / (2. * (np.cos(m * np.pi / s) + np.cos(n * np.pi / s) - 2.)) 
#phi_fft = kappa_fft / (np.sin(m * np.pi / s)**2 + np.sin(n * np.pi / s)**2) 
phi = np.fft.ifft2(np.fft.ifftshift(phi_fft)).real

# k = 2pi / lambda
# 

kappa_check = np.fft.ifft2(kappa_fft).real

pylab.figure()
#pylab.pcolor(kappa_fft.real, cmap='binary')
#pylab.pcolor(phi, cmap='binary')
pylab.contourf(n, m, phi, 100, cmap='binary')
pylab.colorbar()
#pylab.xlim(0., s)
#pylab.ylim(0., s)
pylab.xlim(-s / 2, s / 2)
pylab.ylim(-s / 2, s / 2)

rho = 0.5 * (np.gradient(np.gradient(phi, axis=0), axis=0) + np.gradient(np.gradient(phi, axis=1), axis=1))
#rho += (np.sum(kappa) - np.sum(rho)) / s**2

pylab.figure()
#pylab.pcolor(np.gradient(phi, axis=0), cmap='binary')
#pylab.contourf(n, m, np.gradient(phi, axis=0), 100, cmap='binary')
pylab.contourf(n, m, rho, 100, cmap='binary')
pylab.colorbar()
#pylab.xlim(0., s)
#pylab.ylim(0., s)
pylab.xlim(-s / 2, s / 2)
pylab.ylim(-s / 2, s / 2)

# Now do the ray-tracing

n_trace = 10000000
x = np.random.uniform(-s / 2, s / 2, n_trace)
y = np.random.uniform(-s / 2, s / 2, n_trace)

#x = np.random.uniform(-10., 10., n_trace)
#y = np.random.uniform(-10., 10., n_trace)

grad_x = np.gradient(phi, axis=0)
grad_y = np.gradient(phi, axis=1)

index_x = np.digitize(x, bins=np.arange(-s / 2, s / 2)) - 1
index_y = np.digitize(y, bins=np.arange(-s / 2, s / 2)) - 1

#f_x = scipy.interpolate.interp2d(n, m, np.gradient(phi, axis=0), kind='linear')
#f_y = scipy.interpolate.interp2d(n, m, np.gradient(phi, axis=1), kind='linear')

x_source = x + grad_x[index_x, index_y]
y_source = y + grad_y[index_x, index_y]

#h2 = np.hist2d(x_source, y_source, bins=(np.arange(-s / 2, s / 2, np.arange(-s / 2, s / 2)))[0]
h2 = np.histogram2d(x_source, y_source, bins=np.arange(-s / 2, 1 + s/2))[0]
h2 /= n_trace / (h2.shape[0] * h2.shape[1]) 

pylab.figure()
#pylab.scatter(x, y, c=grad_y[index_x, index_y], edgecolor='none', marker='o')
pylab.contourf(n, m, h2, 20)
pylab.colorbar()
