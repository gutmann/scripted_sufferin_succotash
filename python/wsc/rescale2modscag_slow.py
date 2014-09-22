#!/usr/bin/env python
# (c) Ethan Gutmann gutmann@ucar.edu
from __future__ import print_function
import numpy as np
import mygis

from bunch import Bunch

from wsc import modscag, snodas, wrf, nlcd, regrid
from stat_down import regrid_hi2low as agg

def load_hi_res_dem():
    """docstring for load_hi_res_dem"""
    dem_file="/d2/gutmann/wsc/dem/DEM_CO_NHDPlus_1_arc_seconds.nc"
    x0=2500;x1=4000;y0=12000;y1=10500
    # x0=0;x1=None;y0=-1;y1=0
    demf=mygis.read_nc(dem_file,"elev_m",returnNCvar=True)
    dem=demf.data[y0:y1:-1,x0:x1]
    demf.ncfile.close
    latf=mygis.read_nc(dem_file,"lat",returnNCvar=True)
    lat=latf.data[y0:y1:-1]
    latf.ncfile.close()
    lonf=mygis.read_nc(dem_file,"lon",returnNCvar=True)
    lon=lonf.data[x0:x1]
    lonf.ncfile.close()
    # dem=mygis.read_nc(dem_file,"elev_m").data[::-1,:]
    # lat=mygis.read_nc(dem_file,"lat").data[::-1]
    # lon=mygis.read_nc(dem_file,"lon").data
    lon,lat=np.meshgrid(lon,lat)
    return Bunch(data=dem,lat=lat,lon=lon)


def main():
    """docstring for main"""
    dem=load_hi_res_dem()
    nx=dem.data.shape[1]
    ny=dem.data.shape[0]
    
    for i in range(ny):
        for j in range(nx):
            curlat=dem.lat[i,j]
            curlon=dem.lon[i,j]
            
    


if __name__ == '__main__':
    main()