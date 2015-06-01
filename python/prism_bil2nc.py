#!/usr/bin/env python

import mygis
import numpy as np

import glob
from bunch import Bunch

def read_geo_header(filename):
    """docstring for read_geo"""
    info=Bunch()
    with open(filename,"r") as f:
        for l in f:
            keyname=l.split()[0]
            value=l.split()[1]
            info[keyname]=value
    
    top_lat=float(info["ULYMAP"])
    left_lon=float(info["ULXMAP"])
    dx=float(info["XDIM"])
    dy=float(info["YDIM"])
    nx=int(info["NCOLS"])
    ny=int(info["NROWS"])
    
    info.ny=ny
    info.nx=nx
    info.dy=dy
    info.dx=dx
    info.top_lat=top_lat
    info.left_lon=left_lon
    
    south_lat=top_lat-(ny-1)*dy #-1 because top_lat is the end point
    north_lat=top_lat+dy #+dy because end point has to be past the last value you want
    west_lon=left_lon
    east_lon=left_lon+nx*dx
    
    info.lat=np.arange(south_lat, north_lat, dy)
    info.lon=np.arange(west_lon,  east_lon,  dx)
    
    if (info.lat.shape[0]!=ny) or (info.lon.shape[0]!=nx):
        raise ValueError("nx,ny don't match geo specs")
    
    return info

def read_band_interlaced_by_line(filename,info):
    """docstring for read_band_interlaced_by_line"""
    data=np.fromfile(filename,dtype=np.float32) # read in the raw data
    data=data.reshape((info.ny,info.nx))     # change the shape to a 2d array
    data=data[::-1]                          # reverse the order on the y axis
    return data

def main():
    """docstring for main"""
    files = glob.glob("*_bil.bil")
    info = read_geo_header(glob.glob("*.hdr")[0])
    
    for f in files:
        print("-"*30)
        print("Reading:"+f)
        data = read_band_interlaced_by_line(f,info)
        outputfile=f.replace("_bil.bil",".nc")
        print("Writing:"+outputfile)
        atts=Bunch(units="mm")
        atts["_FillValue"]=info["NODATA"]
        mygis.write(outputfile, data,  varname="pr", dims=("lat","lon"), attributes=atts,
                extravars=[Bunch(data=info.lat,name="lat",dtype="f",dims=("lat",),attributes=Bunch(units="degrees")),
                           Bunch(data=info.lon,name="lon",dtype="f",dims=("lon",),attributes=Bunch(units="degrees"))])
    
    
if __name__ == '__main__':
    main()