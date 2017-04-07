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

import instcat


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
                            help='distance to satellite (kpc)')

        group = parser.add_argument_group('Kernel')
        group.add_argument('--kernel',type=str,default='Gaussian',
                           help='kernel type')
        group.add_argument('--ra',type=float,default=54.0,
                           help='centroid right acension (deg)')
        group.add_argument('--dec',type=float,default=-54.0,
                           help='centroid declination (deg)')
        group.add_argument('--position_angle',type=float,default=0.0,
                           help='position angle east-of-north (deg)')
        # ignore Position Angle for now, assuming the stream is east-west oriented.


        # this code will not using the following parameters:
        # stellar_mass
        # absolute_magnitude
        # extension
        # ellipticity
        # half_light_radius

        # Extra terms for streams.
        group = parser.add_argument_group('Additional parameters for streams')
        group.add_argument('--surface_brightness',type=float,default=30,
                           help='average surface brightness within the width of stream (mag/arcsec^2)')
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--angular_width',type=float,default=0.3,
                           help='angular width (FWHM) of stream (deg)')
        egroup.add_argument('--width',type=float,default=None,
                           help='physical width (FWHM) of stream (pc)')
        egroup = group.add_mutually_exclusive_group()
        egroup.add_argument('--angular_length',type=float,default=3.,
                           help='angular length of stream (deg)')
        egroup.add_argument('--length',type=float,default=None,
                           help='physical length of stream (pc)')
        # assuming Gaussian across stream.
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

    print 'distance_modulus', distance_modulus

    isochrone=Dwarf.createIsochrone(name=args.isochrone, age=args.age,
                                    metallicity=args.metallicity,
                                    distance_modulus=distance_modulus)
    dwarf.set_isochrone(isochrone)

    if args.width == None:
        angular_width = args.angular_width
    else:
        # Convert physical width to angular width in degrees
        angular_width = args.width/(args.distance*1000) * 180./np.pi

    if args.length == None:
        angular_length = args.angular_length
    else:
        # Convert physical length to angular width in degrees
        angular_length = args.length/(args.distance*1000) * 180./np.pi

    print 'angular_width:', angular_width, 'deg'
    print 'angular_length:', angular_length, 'deg'

    angular_radius = angular_width / 2.35 # in degree
    print 'angular_radius:', angular_radius, 'deg'

    ndwarf  = int(angular_length/angular_radius)
    print 'going to simulate', ndwarf, 'dwarfs'

    #convert from surface brightness to absolute magnitude
    #this is a simple estimation and need to be modified in the future
    area = angular_width * 3600 * angular_length * 3600
    apparent_magnitude = args.surface_brightness - 2.5 * np.log10(area)
    absolute_magnitude = apparent_magnitude - distance_modulus
    print 'stream extension:', area, 'in arcsec^2'
    print 'surface brightness', args.surface_brightness
    print 'total apparent magnitude', apparent_magnitude
    print 'total absolute magnitude', absolute_magnitude

    absolute_magnitude_single = absolute_magnitude + 2.5 * np.log10(ndwarf)
    print 'single absolute magnitude', absolute_magnitude_single

    #convert from absolute magnitude to stellar mass / richness
    from scipy.interpolate import UnivariateSpline
    rich = np.logspace(2., 9., 1000)
    mag = isochrone.absolute_magnitude(rich)
    rich = rich[np.argsort(mag)]
    mag = np.sort(mag)
    mag_to_rich = UnivariateSpline(mag, rich, s=0.)
    dwarf.richness = mag_to_rich(absolute_magnitude_single)

    #generate array for the center of ndwarfs
    ra_arr = args.ra + np.arange(-(ndwarf-1.)/2, (ndwarf-1.)/2+1., 1.) * angular_radius / np.cos(np.deg2rad(args.dec))
    dec_arr = args.dec + np.zeros(ndwarf)


    writer = instcat.InstCatWriter()

    data_all = [] # initial the big dataset for writing
    k = ndwarf #used for reassign the object ID, starting with ndwarf because 0-ndwarf will be for unresolved background

    # Write output
    if args.outfile:
        outfile = args.outfile
        if args.nsims > 1:
            base, ext = os.path.splitext(outfile)
            suffix = '_{:0{width}d}'.format(i + 1, width=len(str(args.nsims)))
            outfile = base + suffix + ext
        if os.path.exists(outfile): os.remove(outfile)
        logging.info("Writing %s..." % outfile)
        out = open(outfile, 'w', 1)
    else:
        out = sys.stdout

    kernel = Dwarf.createKernel(name=args.kernel, extension=angular_radius,
                                ellipticity=0,
                                position_angle=0,
                                lon=args.ra, lat=args.dec)
    dwarf.set_kernel(kernel)

    writer.write_toppart(out, dwarf)

    print 'center ra, dec:', dwarf.lon, dwarf.lat

    for j in range(args.nsims):
        for i in range(ndwarf):
            kernel=Dwarf.createKernel(name=args.kernel,extension=angular_radius,
                                      ellipticity=0,
                                      position_angle=0,
                                      lon=ra_arr[i],lat=dec_arr[i])
            dwarf.set_kernel(kernel)


            # Set the kinematic properties
            if args.rs is not None: args.rvmax = 2.163*args.rs
            if args.rhos is not None: raise Exception('Not implemented')
            kinematics=Dwarf.createKinematics(name=args.kinematics,
                                              vdisp=args.vdisp, vmean=args.vmean,
                                              vmax=args.vmax, rvmax=args.rvmax)
            dwarf.set_kinematics(kinematics)
            logging.debug(str(dwarf))


            # Run the simulation
            logging.info("Simulating galaxy %i..."%i)
            data = Simulator.simulate(dwarf)
            data['ID'] = data['ID'] + k
            k = data['ID'][-1]
            writer.write_middlepart(out,dwarf,data,i)
            data_all.append(data)

        print 'each dwarf has', len(data), 'stars'
        print 'the stream has', len(data_all)*len(data), 'stars'

        data_all = np.concatenate(data_all)


        writer.write_endpart(out,dwarf,data_all)
        out.flush()
