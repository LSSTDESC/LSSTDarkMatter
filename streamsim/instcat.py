#!/usr/bin/env python
"""
Convert a simulated set of stars into an instcat file for imsim.
"""
__author__ = "Alex Drlica-Wagner"

import os
import copy
from collections import OrderedDict as odict
import numpy as np

FILTER = 'g'
MAG = 'MAG_%s'%FILTER.upper()

MAGLIMS =odict([('u',23),('g',25),('r',25),('i',24.5),('z',24),('y',23)])
FILTERS = odict([('u',0),('g',1),('r',2),('i',3),('z',4),('y',5)])

DEFAULTS = odict([
    ('rightascension', 31.1133844),
    ('declination', -10.0970060),
    ('mjd', 59797.2854090),
    ('altitude', 43.6990272),
    ('azimuth', 73.7707957),
    ('filter', 2),
    ('rotskypos', 69.0922930),
    ('rottelpos', 0.000),
    ('dist2moon', 145.1095257),
    ('moonalt', -11.1383568),
    ('moondec', -18.7702120),
    ('moonphase', 59.6288830),
    ('moonra', 230.9832941),
    ('nsnap', 2),
    ('obshistid', 1),
    ('seed', 1),
    ('seeing', 0.7613760),
    ('rawSeeing', 0.7613760),
    ('FWHMgeom', 1.210734),
    ('FWHMeff', 1.409653),
    ('sunalt', -59.1098785),
    ('vistime', 33.0000000),
    ])

POINT = "object %(ID)i %(RA)12.6f %(DEC)12.6f %({})12.6f starSED/phoSimMLT/lte031-4.5-1.0a+0.4.BT-Settl.spec.gz 0 0 0 0 0 0 point none CCM 0.0 3.1".format(MAG)

SERSIC = "object %(ID)i %(RA)12.6f %(DEC)12.6f %({})12.6f starSED/phoSimMLT/lte031-4.5-1.0a+0.4.BT-Settl.spec.gz 0 0 0 0 0 0 sersic2d %(EXT)12.6f %(EXT)12.6f 0 0.5 none CCM 0.0 3.1".format(MAG)

def set_defaults(kwargs,defaults):
    for k,v in defaults.items():
        kwargs.setdefault(k,v)
    return kwargs

def read_instcat(filename,**kwargs):
    """ Read an instcat file into a numpy.recarray """
    defaults =dict(skip_header=23,usecols=[1,2,3,4],
                   dtype=[('ID',int),('RA',float),('DEC',float),('MAG',float)])
    set_defaults(kwargs,defaults)
    return np.genfromtxt(filename,**kwargs)


class InstCatWriter(object):

    def write_toppart(self, filename, dwarf, force=True):
        if isinstance(filename,basestring):
            if os.path.exists(filename) and not force:
                msg = "File %s already exists"%filename
                raise IOError(msg)
            out = open(filename,'w')
        elif isinstance(filename,file):
            out = filename
        else:
            msg ="Unrecongnized file object"
            raise IOError(msg)

        out.write(self.comment_string())
        out.write(self.header_string(dwarf))
        return out


    def write_middlepart(self, filename, dwarf, data, idx, force=True):
        if isinstance(filename,basestring):
            if os.path.exists(filename) and not force:
                msg = "File %s already exists"%filename
                raise IOError(msg)
            out = open(filename,'a')
        elif isinstance(filename,file):
            out = filename
        else:
            msg ="Unrecongnized file object"
            raise IOError(msg)
        out.write(self.unresolved_string(dwarf,data, idx))

        return out

    def write_endpart(self, filename, dwarf, data, force=True):
        if isinstance(filename,basestring):
            if os.path.exists(filename) and not force:
                msg = "File %s already exists"%filename
                raise IOError(msg)
            out = open(filename,'a')
        elif isinstance(filename,file):
            out = filename
        else:
            msg ="Unrecongnized file object"
            raise IOError(msg)
        out.write(self.resolved_string(dwarf,data))
        out.write('\n')
        return out

    def comment_string(self,comment=None):
        if comment is None: comment = "instcat file written by dsphsim"
        return "# %s \n"%comment

    def header_string(self, dwarf,**kwargs):
        header = self.parse_dwarf(dwarf,**kwargs)
        return "\n".join(["%s %s"%(k,v) for k,v in header.items()])+'\n'

    def parse_dwarf(self, dwarf, **kwargs):
        header = set_defaults(odict(),DEFAULTS)
        header.update(kwargs)
        header['rightascension'] = dwarf.lon
        header['declination'] = dwarf.lat
        header['filter'] = FILTERS[FILTER]
        return header

    def resolved_string(self, dwarf, data, **kwargs):
        sel = (data[MAG] <= MAGLIMS[FILTER])
        return '\n'.join(np.char.mod([POINT],data[sel]))

    def unresolved_string(self, dwarf, data, idx, **kwargs):
        sel = (data[MAG] > MAGLIMS[FILTER])
        params = dict(ID=idx,
                      RA=dwarf.lon,
                      DEC=dwarf.lat,
                      EXT=dwarf.extension*3600 / 1.68, # arcsec
        )
        params[MAG] = self.cumulative_magnitude(data[sel])
        return SERSIC%(params) + '\n'

    def cumulative_magnitude(self,data):
        return -2.5 * np.log10(np.sum(10**(data[MAG]/-2.5)))

