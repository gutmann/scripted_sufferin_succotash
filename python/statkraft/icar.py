from __future__ import absolute_import, print_function, division

import numpy as np

import mygis
from bunch import Bunch

def geo(filename):
    lats=mygis.read_nc(filename,"lat").data
    lons=mygis.read_nc(filename,"lon").data

    return Bunch(lon=lons,lat=lats)
