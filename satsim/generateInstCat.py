#!/usr/bin/env python
"""
Simple tool to simulate stellar populations.
"""

import os,sys
import numpy as np
import logging

import scipy.stats as stats

import dsphsim
from dsphsim.dwarf import Dwarf
from dsphsim.instruments import factory as instrumentFactory
from dsphsim.tactician import factory as tacticianFactory
from dsphsim.velocity import PhysicalVelocity

class Simulator(object):

    @staticmethod
    def simulate(dwarf,**kwargs):
        """ Simulate dwarf galaxy """

        # Set the second band to 'i' (matches CaT lines)
        dwarf.band_1 = 'g'; dwarf.band_2 = 'r'
        mag_1,mag_2,ra,dec,velocity = dwarf.simulate()
        angsep = dwarf.kernel.angsep(ra,dec)
        rproj = dwarf.distance * np.tan(np.radians(angsep))
        idx = np.arange(1,len(mag_1)+1)

        # Do we also want to save vsyserr as VSYSERR?
        names = ['ID','RA','DEC',
                 'MAG_%s'%dwarf.band_1.upper(),'MAG_%s'%dwarf.band_2.upper(),
                 'ANGSEP','RPROJ']
        data = [idx, ra, dec,
                mag_1, mag_2,
                angsep, rproj]
        return np.rec.fromarrays(data,names=names)

    @classmethod
    def parser(cls):
        import argparse
        description = "Simulate the observable properties of a dwarf galaxy."
        formatter = argparse.ArgumentDefaultsHelpFormatter
        parser = argparse.ArgumentParser(description=description,
                                         formatter_class=formatter)
        parser.add_argument('outfile',nargs='?',
                            help="optional output file")
        parser.add_argument('--seed',type=int,default=None,
                            help="random seed")
        parser.add_argument('-v','--verbose',action='store_true',
                            help="verbose output")
        parser.add_argument('-n','--nsims',type=int,default=1,
                            help="number of simulations")

        group = parser.add_argument_group('Physical')
        # Stellar mass and absolute magnitude are mutually exclusive
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--stellar_mass',type=float,default=2000.,
                            help='stellar mass of satellite (Msun)')
        egroup.add_argument('--absolute_magnitude',type=float,default=None,
                            help='absolute magnitude of satellite')


        group = parser.add_argument_group('Kinematic')
        group.add_argument('--kinematics',type=str,default='Gaussian',
                           help='kinematic distribution function')
        group.add_argument('--vmean',type=float,default=60.,
                            help='mean systemic velocity (km/s)')
        # should be mutually exclusive with vmax and rs
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--vdisp',type=float,default=3.3,
                            help='gaussian velocity dispersion (km/s)')
        egroup.add_argument('--vmax',type=float,default=10.0,
                            help='maximum circular velocity (km/s)')
        egroup.add_argument('--rhos',type=float,default=None,
                            help='maximum circular velocity (Msun/pc^3)')
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--rvmax',type=float,default=0.4,
                           help='radius of max circular velocity (kpc)')
        # ADW: it would be nice to remove this...
        egroup.add_argument('--rs',type=float,default=None,
                           help='scale radius for NFW halo (kpc)')

        group = parser.add_argument_group('Isochrone')
        group.add_argument('--isochrone',type=str,default='Bressan2012',
                            help='isochrone type')
        group.add_argument('--age',type=float,default=12.0,
                           help='age of stellar population (Gyr)')
        group.add_argument('--metallicity',type=float,default=2e-4,
                           help='metallicity of stellar population')
        # Distance modulus and distance are mutually exclusive
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--distance_modulus',type=float,default=17.5,
                            help='distance modulus')
        egroup.add_argument('--distance',type=float,default=None,
                            help='distaance to satellite (kpc)')

        group = parser.add_argument_group('Kernel')
        group.add_argument('--kernel',type=str,default='RadialExponential',
                           help='kernel type')
        group.add_argument('--ra',type=float,default=54.0,
                           help='centroid right acension (deg)')
        group.add_argument('--dec',type=float,default=-54.0,
                           help='centroid declination (deg)')
        group.add_argument('--extension',type=float,default=0.1,
                           help='projected half-light radius (deg)')
        group.add_argument('--ellipticity',type=float,default=0.0,
                           help='spatial extension (deg)')
        # Physical and angular extension are mutually exclusive
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--position_angle',type=float,default=0.0,
                           help='position angle east-of-north (deg)')
        egroup.add_argument('--half_light_radius',type=float,default=None,
                            help='half-light radius (pc)')

        return parser

if __name__ == "__main__":
    parser = Simulator.parser()
    args = parser.parse_args()
    kwargs = vars(args)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.seed is not None:
        np.random.seed(args.seed)

    dwarf = Dwarf()
    
    if args.distance == None:
        distance_modulus = args.distance_modulus
    else:
        # Calculate distance modulus from distance in kpc
        distance_modulus = 5. * (np.log10(args.distance * 1.e3) - 1.)

    isochrone=Dwarf.createIsochrone(name=args.isochrone, age=args.age,
                                    metallicity=args.metallicity,
                                    distance_modulus=distance_modulus)
    dwarf.set_isochrone(isochrone)

    if args.half_light_radius == None:
        extension = args.extension
    else:
        # Convert physical half-light radius to extension in degrees
        extension = args.half_light_radius/(args.distance*1000) * 180./np.pi

    kernel=Dwarf.createKernel(name=args.kernel,extension=extension,
                              ellipticity=args.ellipticity,
                              position_angle=args.position_angle,
                              lon=args.ra,lat=args.dec)
    dwarf.set_kernel(kernel)
    
    if args.absolute_magnitude == None:
        stellar_mass = args.stellar_mass
        dwarf.richness = stellar_mass/dwarf.isochrone.stellar_mass()
    else:         
        from scipy.interpolate import UnivariateSpline
        rich = np.logspace(2., 9., 1000)
        mag = isochrone.absolute_magnitude(rich)
        rich = rich[np.argsort(mag)]
        mag = np.sort(mag)
        mag_to_rich = UnivariateSpline(mag, rich, s=0.)
        dwarf.richness = mag_to_rich(args.absolute_magnitude)
        # print 'Recalculated abs mag = %.2f' %isochrone.absolute_magnitude(dwarf.richness)
        # print 'Input abs mag = %.2f' %args.absolute_magnitude

    # Set the kinematic properties
    if args.rs is not None: args.rvmax = 2.163*args.rs
    if args.rhos is not None: raise Exception('Not implemented')
    kinematics=Dwarf.createKinematics(name=args.kinematics,
                                      vdisp=args.vdisp, vmean=args.vmean,
                                      vmax=args.vmax, rvmax=args.rvmax)
    dwarf.set_kinematics(kinematics)
    logging.debug(str(dwarf))

    import instcat
    writer = instcat.InstCatWriter()

    for i in range(args.nsims):
        # Run the simulation
        logging.info("Simulating galaxy %i..."%i)
        data = Simulator.simulate(dwarf)

        # Write output
        if args.outfile:
            outfile = args.outfile
            base,ext = os.path.splitext(outfile)
            if args.nsims > 1:
                suffix = '_{:0{width}d}'.format(i+1,width=len(str(args.nsims)))
                outfile_g = base + '_g' + suffix + ext
                outfile_r = base + '_r' + suffix + ext
            else:
                outfile_g = base + '_g' + ext
                outfile_r = base + '_r' + ext
            print ('outfile_g = %s' %outfile_g)
            print ('outfile_r = %s' %outfile_r)
            if os.path.exists(outfile_g): os.remove(outfile_g)
            if os.path.exists(outfile_r): os.remove(outfile_r)
            logging.info("Writing %s and %s..."%(outfile_g, outfile_r))
            out_g = open(outfile_g,'w',1)
            out_r = open(outfile_r,'w',1)
        else:
            out_g = sys.stdout
            out_r = sys.stdout

        writer.write_g(out_g,dwarf,data)
        writer.write_r(out_r,dwarf,data)
        out_g.flush()
        out_r.flush()
