#!/usr/bin/env python
from __future__ import print_function
import glob
import sys
import re

import numpy as np

# NOTE, you need these libraries in your PYTHONPATH
import shapefile
from point_in_poly import points_in_poly
from stat_down import myio as io

def load_shapes(filename):
    """Read in all Shape Records from a given shapefile"""
	shp=shapefile.Reader(filename)
	shprecs=shp.shapeRecords()
	return shprecs

def load_geo(filename,latvar=None,lonvar=None):
    """Load Geographic (lat/lon coords) from a netCDF file
    
    Output will be a (lat,lon) tuple of 2d grids with lon=[-180 - 180]
    """
    lat,lon=None,None
    # if lat/lon var names were specified just read those data
    if latvar:
        try:
            lat=io.read_nc(filename,latvar).data
            lon=io.read_nc(filename,lonvar).data
        except exception as e:
            print(e)
    # otherwise you have to search a list of 
    # variable names to test
    latnames=["lat","latitude","XLAT"]
    lonnames=["lon","longitude","XLONG"]
    # loop over possible variable names
    i=0
    while lat==None:
        try:
            lat=io.read_nc(filename,latnames[i]).data
            lon=io.read_nc(filename,lonnames[i]).data
        except:
            lat=None
            lon=None
        i+=1
    # if this grid is 0-360 convert it to -180 - 180
    # note this actually converts to the opposite sign convention as normal...
    if lon.max()>180:
        lon-=360
    # if grids are 1D make them 2D
    if len(lat.shape)==1:
        lon,lat=np.meshgrid(lon,lat)
    
    return lat,lon


def write_file(filename,data,lat,lon):
    """write the result as a 2D netCDF grid"""
    io.write(filename,data.astype("d"),dtype='d',lat=lat,lon=lon)


def mark_poly(data,poly,lat,lon,code_position):
    """find all data points that are in a given polygon and mark them with the current HUC
    data = grid to mark
    poly = a shapefile polygon dataset
    lat,lon = 2d grids of latitude and longitude
    code_position = index into the polygon record info that stores the HUC code
    """
    polynumber=poly.record[code_position]
    # the bounding box to VASTLY speed up calculations
    bbox=poly.shape.bbox
    # we need to search each part independantly
    parts=poly.shape.parts
    points=poly.shape.points
    # first find the points that are in the bounding box
    inbox=np.where((lon>bbox[0]) & (lon<=bbox[2]) & (lat>=bbox[1]) & (lat<=bbox[3]))
    for i in range(len(parts)):
        # subset the polygon points for the current part
        if i==len(parts)-1:
            xypoints=np.array(points[parts[i]:])
        else:
            xypoints=np.array(points[parts[i]:parts[i+1]])
        # find the points that fall inside the current polygon part
        internalpoints=points_in_poly(lon[inbox],lat[inbox],xypoints)
        # mark those points with this HUC code
        data[inbox[0][internalpoints],inbox[1][internalpoints]]=polynumber

def main(geo_file,poly_file,outputfilename):
    print("Loading Shapefile")
    shapes=load_shapes(poly_file)
    print("Loading Geographic data")
    lat,lon=load_geo(geo_file)
    
    # variable to hold the outputdata
    outputdata=np.zeros(lat.shape,dtype=np.long)
    code_position=11
    # loop over all shapes in the shapefile
    for i in range(len(shapes)):
        # print a progress report
        print("Progress={0}/{1}       ".format(i, len(shapes)),end="\r")
        sys.stdout.flush()
        # here is where we do the actual work
        mark_poly(outputdata,shapes[i],lat,lon,code_position)
    print("\nFinished")
    write_file(outputfilename,outputdata,lat,lon)
    
if __name__ == '__main__':
    # stupidly simple commandline processing,
    # use argparse if we need more options in the future. 
    try:
        main(sys.argv[1],sys.argv[2],sys.argv[3])
    except Exception as e:
        print(e)
        print("USAGE:")
        print("""  hucmask.py <rasterfilename.nc> <shapefilename.shp> <outputfilename>
        
        rasterfile must have lat/lon, latitude/longitude, or XLAT/XLONG variables
        shapefile needs an associated dbf file in the same directory
        """)
