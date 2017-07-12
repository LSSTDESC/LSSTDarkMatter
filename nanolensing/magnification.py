import numpy as np
import pylab

pylab.ion()

kappa_array = np.linspace(0., 2., 1000)
shear_array = np.linspace(0., 1., 1000)

kappa_mesh, shear_mesh = np.meshgrid(kappa_array, shear_array)

delta_kappa_mesh = np.tile(1.e-6, kappa_mesh.shape)

mu_mesh = 1. / np.fabs((1 - kappa_mesh)**2 - shear_mesh**2)

delta_mu_mesh = 1. / np.fabs((1 - (kappa_mesh + delta_kappa_mesh))**2 - shear_mesh**2)

delta_mesh = (delta_mu_mesh - mu_mesh) / delta_kappa_mesh

pylab.figure()
pylab.contourf(kappa_mesh, shear_mesh, mu_mesh, levels=np.linspace(1., 100., 20), cmap='Blues_r', extend='max')
#pylab.contourf(kappa_mesh, shear_mesh, delta_mesh, levels=np.linspace(-100., 100., 50), cmap='jet')
colorbar = pylab.colorbar()
colorbar.set_label('Magnification')
#pylab.scatter(0.392, 0.642, c='white', s=200)
pylab.xlabel('Convergence')
pylab.ylabel('Shear')
pylab.xlim(0., 2.)
pylab.ylim(0., 1.)
pylab.savefig('magnification.pdf')

#x_mesh, y_mesh = np.meshgrid(np.arange(10), np.arange(0, 100, 10))
#z_mesh = x_mesh + y_mesh
#pylab.figure()
#pylab.contourf(x_mesh, y_mesh, z_mesh, levels=np.arange(100))
#pylab.colorbar()
