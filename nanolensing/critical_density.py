import numpy as np
import pylab

import cosmology
import cosmology_constants

pylab.ion()

#####

z_array = np.linspace(0., 10., 10000)

pylab.figure()
pylab.plot(z_array, cosmology.D_A(z_array) / cosmology.D_H)

#####

z_array = np.linspace(0., 5., 100)
z_1, z_2 = np.meshgrid(z_array, z_array)

d_a_12 = np.clip(cosmology.D_A(z_1, z_2) / cosmology.D_H, 0., np.inf)

pylab.figure()
pylab.contourf(z_1, z_2, d_a_12, 100, cmap='jet')
pylab.colorbar()
pylab.xlabel('z_lens')
pylab.ylabel('z_source')

#####

distance_term = cosmology.D_A(z_2) / (np.clip(cosmology.D_A(z_1, z_2), 0., np.inf) * cosmology.D_A(z_1))
critical_density = cosmology_constants.C_LIGHT**2 * distance_term / (4. * np.pi * cosmology_constants.G_NEWTON) # g cm^-2
critical_density *= (cosmology_constants.M_SOLAR_TO_G)**(-1) * cosmology_constants.PC_TO_CM**2

print np.min(critical_density), np.max(critical_density) 


pylab.figure()
#pylab.pcolor(critical_density, vmin=0., vmax=10000)
pylab.contourf(z_1, z_2, np.clip(critical_density, 0., 10000), 100, cmap='jet')
#pylab.contourf(z_1, z_2, cosmology.D_A(z_1) / cosmology.D_H, 100, cmap='jet')
#pylab.contourf(z_1, z_2, c=critical_density, cmap='jet', edgecolor='none')
pylab.colorbar()
pylab.xlabel('z_lens')
pylab.ylabel('z_source')

#z_1 = 0.3
#z_2 = 0.31
#cosmology.D_A(z_2) / (cosmology.D_A(z_1, z_2) * cosmology.D_A(z_1))

"""
# Check against https://arxiv.org/pdf/1103.0353.pdf
z_1 = 0.96
z_2 = 2.64
distance_term = cosmology.D_A(z_2) / (np.clip(cosmology.D_A(z_1, z_2), 0., np.inf) * cosmology.D_A(z_1))
critical_density = cosmology_constants.C_LIGHT**2 * distance_term / (4. * np.pi * cosmology_constants.G_NEWTON) # g cm^-2
critical_density *= (cosmology_constants.M_SOLAR_TO_G)**(-1) * cosmology_constants.PC_TO_CM**2 # M_sol pc^-2
print critical_density
"""
