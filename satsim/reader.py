#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np
from writer import set_defaults

def read_instcat(filename,**kwargs):
    defaults =dict(skip_header=23,usecols=[1,2,3,4],
                   dtype=[('ID',int),('RA',float),('DEC',float),('MAG',float)])
    set_defaults(kwargs,defaults)

    return np.genfromtxt(filename,**kwargs)

