#!/usr/bin/env python
from __future__ import print_function
import gc,sys
import numpy as np

import mygis
from bunch import Bunch

import snodas

# from flushprint import Flushfile
# sys.stdout=Flushfile(sys.stdout)

dem_file="/d5/gutmann/lidar/grass/snow-off-dem.tif"
nlcd_file="/d2/gutmann/nldc/NLCD2006_landcover.tif"
snodas_file="/d2/gutmann/wsc/snodas/SWE_Daily0600UTC_WesternUS_2010.dat"
nlcd_subset=[25000,75000,30000,63000]
nlcd_subset=[40000,52000,52000,58500]
snodas_subset=[0,None,1200,None]

# nlcd_subset=[42000,47000,54000,58000]
# nlcd_subset=[44000,46000,55000,56000]
# snodas_subset=[850,900,1980,2030]

def geo2lat_lon(data):
    """calculate xy points in data projection then convert to lat/lon"""
    left=data.geo[0]
    top=data.geo[3]
    dx=data.geo[1]
    dy=data.geo[5]
    
    nx=data.data.shape[1]
    ny=data.data.shape[0]
    
    x=np.arange(left,left+dx*nx,dx)[nlcd_subset[2]:nlcd_subset[3]].astype(np.int32)
    y=np.arange(top, top +dy*ny,dy)[nlcd_subset[0]:nlcd_subset[1]].astype(np.int32)

    return Bunch(x=x,y=y)
    
    # nx=nlcd_subset[3]-nlcd_subset[2]
    # ny=nlcd_subset[1]-nlcd_subset[0]
    # x,y=np.meshgrid(x,y)
    # x=x.reshape((nx*ny))
    # y=y.reshape((nx*ny))
    # 
    # points=np.concatenate([x[:,np.newaxis],y[:,np.newaxis]],axis=1)
    # 
    # print("Converting Points from albers to lat/lon")
    # lat,lon=mygis.proj2ll(points=points,proj=data.proj)
    # lat=lat.reshape((ny,nx))
    # lon=lon.reshape((ny,nx))
    # return Bunch(lat=lat,lon=lon)

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

def lower_resolution(forest_cover,fgeo,ogeo,proj):
    """convert forest cover to a lower resolution map"""
    if len(ogeo.lat.shape)>1:
        olon,olat=(ogeo.lon[1,:],ogeo.lat[:,1])
    else:
        olon,olat=(ogeo.lon,ogeo.lat)
    
    outputdata=np.zeros((2,olat.shape[0],olon.shape[0]))
    for i in range(forest_cover.shape[0]):
        if (i%100)==0:
            print(i,forest_cover.shape[0], "    "+str(float(i)/forest_cover.shape[0]*100)[:5]+"%", end=" ")
        for j in range(forest_cover.shape[1]):
            xlat,xlon=mygis.proj2ll(x=float(fgeo.x[j]),y=float(fgeo.y[i]),proj=proj)
            xpos=np.argmin(np.abs(olon-xlon))
            ypos=np.argmin(np.abs(olat-xlat))
            outputdata[0,ypos,xpos]+=forest_cover[i,j]
            outputdata[1,ypos,xpos]+=1
            
    forest_fraction=outputdata[0,...]/outputdata[1,...]
    forest_fraction[forest_fraction>=0.5]=1
    forest_fraction[forest_fraction<0.5]=0
    return forest_fraction

def forest_cover(geo_data=None):
    print("Loading NLCD file")
    lc_data=mygis.read_tiff(nlcd_file)
    print("setting up geo data")
    lc_geo=geo2lat_lon(lc_data)
    print("first subset")
    sub_data1=lc_data.data[nlcd_subset[0]:nlcd_subset[1],nlcd_subset[2]:nlcd_subset[3]]
    # del lc_data.data
    # gc.collect()
    print(sub_data1.shape)
    if geo_data==None:
        print("Loading SNODAS data")
        geo_data=snodas.load(None)
        subset=snodas_subset
    # print("Subset 2")
    subdata=sub_data1
    # subdata=subset_a2b(sub_data1,lc_geo,geo_data,subset)
    # del sub_data1
    # gc.collect()
    
    print("Creating high res forest mask")
    lc_forest=np.zeros(subdata.shape)
    lc_forest[np.where((subdata>40)&(subdata<44))]=1
    
    print("Lowering resolution")
    output_data=lower_resolution(lc_forest,lc_geo,geo_data,lc_data.proj)

    print("Writing data file")
    mygis.write("snodas_forest.nc",output_data)
    
if __name__ == '__main__':
    forest_cover()