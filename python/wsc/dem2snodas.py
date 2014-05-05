#!/usr/bin/env python
from __future__ import print_function
import gc,sys
import numpy as np

import mygis
from bunch import Bunch

import snodas

import flushprint
if flushprint.in_ipython():
    pass
else:
    sys.stdout=flushprint.Flushfile(sys.stdout)

dem_file="/d2/gutmann/wsc/dem/DEM_CO_NHDPlus_1_arc_seconds.nc"
# demgeo_file="/d5/gutmann/lidar/grass/snow-on-dem.nc"
snodas_file="/d2/gutmann/wsc/snodas/SWE_Daily0600UTC_WesternUS_2010.dat"


dims=("lat","lon")
FILL_VALUE=np.float32(-9999) # if this doesn't work, (np.zeros(3)-9999).astype("f")[1] should work...
lat_atts=dict(units="degrees_north",axis="Y", long_name="Latitude",standard_name="latitude")
lat_info=Bunch(data=None,name="lat",dims=("lat",),dtype="f",attributes=lat_atts)

lon_atts=dict(units="degrees_east",axis="X", long_name="Longitude",standard_name="longitude")
lon_info=Bunch(data=None,name="lon",dims=("lon",),dtype="f",attributes=lon_atts)

data_atts=dict(units="m",_FillValue=FILL_VALUE,missing_value=FILL_VALUE)
data_info=Bunch(data=None,name="dem",dims=dims,dtype="f",attributes=data_atts)



def geo2lat_lon(data):
    """calculate xy points in data projection then convert to lat/lon"""
    left=data.geo[0]
    top=data.geo[3]
    dx=data.geo[1]
    dy=data.geo[5]
    
    nx=data.data.shape[1]
    ny=data.data.shape[0]
    
    x=np.arange(left,left+dx*nx,dx)
    y=np.arange(top, top +dy*ny,dy)

    x,y=np.meshgrid(x,y)
    x=x.reshape((nx*ny))
    y=y.reshape((nx*ny))
    
    
    # WARNING: this will break with LARGE numbers of points (Memory error) 
    # e.g. 30m grid over CONUS...
    points=np.concatenate([x[:,np.newaxis],y[:,np.newaxis]],axis=1)
    
    print("Converting Points from data projection to lat/lon")
    lat,lon=mygis.proj2ll(points=points,proj=data.proj)
    lat=lat.reshape((ny,nx))
    lon=lon.reshape((ny,nx))
    return Bunch(lat=lat,lon=lon,x=x,y=y)

def subset_a2b(data,geoin,geoout,subset):
    """subset both data (returned) and geoin to the range of geoout
    OUTPUT/SIDEEFFECTS: 
        geoout is subset to snodas_subset
        geoin  is subset to points within geoout
        subset of data.data is returned
    """
    geoout.lat=geoout.lat[subset[0]:subset[1]]
    geoout.lon=geoout.lon[subset[2]:subset[3]]
    
    minx=np.where(geoin.lon>=geoout.lon.min())[1].min()
    maxx=np.where(geoin.lon<=geoout.lon.max())[1].max()

    # NOTE input lats go from high to low instead of low to high 
    #      the upper left corner of the map is at coordinate (0,0)
    maxy=np.where(geoin.lat>=geoout.lat.min())[0].max()
    miny=np.where(geoin.lat<=geoout.lat.max())[0].min()
    
    geoin.lat=geoin.lat[miny:maxy+1,minx:maxx+1]
    geoin.lon=geoin.lon[miny:maxy+1,minx:maxx+1]
    return data[miny:maxy+1,minx:maxx+1]

def lower_resolution(data,dgeo,ogeo):
    """convert to a lower resolution map by simple averaging"""
    if len(ogeo.lat.shape)>1:
        olon,olat=(ogeo.lon[1,:],ogeo.lat[:,1])
    else:
        olon,olat=(ogeo.lon,ogeo.lat)
    
    outputdata=np.zeros((2,olat.shape[0],olon.shape[0]))
    for i in range(data.shape[0]):
        if (i%100)==0:
            print(str(float(i)/data.shape[0]*100)[:5]+"%", end=",   ")
        for j in range(data.shape[1]):
            xlat,xlon=dgeo.lat[i,j],dgeo.lon[i,j]
            xpos=np.argmin(np.abs(olon-xlon))
            ypos=np.argmin(np.abs(olat-xlat))
            outputdata[0,ypos,xpos]+=data[i,j]
            outputdata[1,ypos,xpos]+=1
            
    return outputdata[0,...]/outputdata[1,...]

def simple_lowres(data,scale=3):
    """docstring for simple_lowres"""
    ny,nx=data.shape
    newny=int(ny/scale)
    newnx=int(nx/scale)
    data=data[:newny*scale,:newnx*scale].reshape((newny,scale,newnx,scale))
    data=data.mean(axis=1).mean(axis=2)

    return data
    
def write(filename,data,geo):
    """write DEM data to a file after subsetting to 'good' data"""
    goodpoints=np.where(data>0)
    ymin=np.min(goodpoints[0])
    ymax=np.max(goodpoints[0])+1
    xmin=np.min(goodpoints[1])
    xmax=np.max(goodpoints[1])+1
    
    output_data=data[ymin:ymax,xmin:xmax]
    output_lat=geo.lat[ymin:ymax]
    output_lon=geo.lon[xmin:xmax]
    
    lat_info.data=output_lat
    lon_info.data=output_lon
    
    data_info.subset="[{},{},{},{}]".format(ymin,ymax,xmin,xmax)
    extra_vars=[lat_info,lon_info]
    mygis.write(filename,output_data,varname=data_info.name,dtype=data_info.dtype,dims=data_info.dims,
               attributes=data_info.attributes,extravars=extra_vars,history="dem2snodas")
    

def rescale():
    print("Loading DEM")
    dem=mygis.read_nc(dem_file,"elev_m").data
    print("setting up geo data")
    dem_geo=mygis.read_geo(dem_file)
    
    dem=simple_lowres(dem)
    dem_geo.lat=simple_lowres(dem_geo.lat)
    dem_geo.lon=simple_lowres(dem_geo.lon)
    
    print("Loading SNODAS coordinates")
    geo_data=snodas.load(None)
    
    print("Lowering resolution")
    output_data=lower_resolution(dem,dem_geo,geo_data)

    print("Writing data file")
    write("snodas_dem.nc",output_data,geo_data)
    
if __name__ == '__main__':
    rescale()