from __future__ import absolute_import, print_function, division

from mpl_toolkits.basemap import pyproj
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

import mygis
from bunch import Bunch

def se_norway_basemap():
    return Basemap(width=500*1000,height=700*1000,resolution="h",lon_0=8.5,lat_0=61,projection="aea")

def vis(data, lat, lon, m=None,
        cmap=plt.cm.gist_rainbow_r, clim=None, title=None,
        cbar=True, lakes=True):
    if m==None:
        m=se_norway_basemap()

    m.pcolormesh(lon,lat,data,latlon=True,cmap=cmap)
    if cbar:
        plt.colorbar()
    if clim!=None:
        plt.clim(clim)
    if title!=None:
        plt.title(title)

    m.drawcoastlines()
    if lakes:
        m.fillcontinents(color=(0,0,0,0.0),lake_color=(0,0,1,1))

    return m
