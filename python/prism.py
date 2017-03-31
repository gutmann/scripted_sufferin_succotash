import numpy as np

from bunch import Bunch

nx=1405
ny=621

def prism_geo():
    ULY = 49.9166666666687
    ULX = -125.0
    dx  = 0.04166666666667
    
    lat = np.arange(ULY - (ny*dx), ULY, dx) + dx
    lon = np.arange(ULX, ULX + (nx*dx), dx)
    
    lon,lat = np.meshgrid(lon,lat)
    
    return Bunch(lon=lon, lat=lat)

def load_bil(filename):
    data = np.fromfile(filename, dtype=np.float32).reshape( (ny,nx) )[::-1,:]
    
    data=np.ma.array(data, mask=(data<0))

    return Bunch(
            data = data,
            geo  = prism_geo() )

def load(filename):
    if filename.split(".")[-1]=="bil":
        return load_bil(filename)
    elif filename.split(".")[-1]=="asc":
        return load_asc(filename)
    else:
        raise IOError
