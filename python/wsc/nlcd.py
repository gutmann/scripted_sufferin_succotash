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

NLCD_FOREST_TYPE=1
NLCD_EXPOSED_TYPE=0
nlcd_file="/d2/gutmann/nldc/NLCD2006_landcover.tif"
snodas_file="/d2/gutmann/wsc/snodas/SWE_Daily0600UTC_WesternUS_2010.dat"
subset_list=[]
subset_list.append([25000,75000,30000,63000]) #large area, takes too long to run
subset_list.append([35000,52000,46000,58500]) # _1 approximate area of CO Rockies
subset_list.append([52000,60000,46000,50000]) # _2 approximate area of CO Rockies
subset_list.append([60000,69000,46000,50000]) # _3 approximate area of CO Rockies
subset_list.append([52000,60000,50000,58000]) # _4 approximate area of CO Rockies
subset_list.append([60000,69000,50000,58000]) # _5 approximate area of CO Rockies
nlcd_subset=[50000,51000,56000,56500] # fast debugging test case
# nlcd_subset=[50000, 55000, 47000, 53000] #san juans
snodas_subset=[0,None,1200,None]

dims=("lat","lon")
FILL_VALUE=np.float32(-9999)
lat_atts=dict(units="degrees_north",axis="Y", long_name="Latitude",standard_name="latitude")
lat_info=Bunch(data=None,name="lat",dims=("lat",),dtype="f",attributes=lat_atts)

lon_atts=dict(units="degrees_east",axis="X", long_name="Longitude",standard_name="longitude")
lon_info=Bunch(data=None,name="lon",dims=("lon",),dtype="f",attributes=lon_atts)

data_atts=dict(units="[]",_FillValue=FILL_VALUE,missing_value=FILL_VALUE,forest="1",exposed="0")
data_info=Bunch(data=None,name="data",dims=dims,dtype="f",attributes=data_atts)


# nlcd_subset=[42000,47000,54000,58000]
# nlcd_subset=[44000,46000,55000,56000]
# snodas_subset=[850,900,1980,2030]

def geo2lat_lon(data,geo_geographic=False,reverse_lat=False):
    """calculate xy points in data projection then convert to lat/lon"""
    left=data.geo[0]
    top=data.geo[3]
    dx=data.geo[1]
    dy=data.geo[5]
    
    nx=data.data.shape[1]
    ny=data.data.shape[0]
    
    x=np.arange(left,left+dx*nx,dx)[nlcd_subset[2]:nlcd_subset[3]].astype(np.int32)
    y=np.arange(top, top +dy*ny,dy)[nlcd_subset[0]:nlcd_subset[1]].astype(np.int32)
    if reverse_lat:
        y=y[::-1]

    if not geo_geographic:
        return Bunch(x=x,y=y)
    
    nx=nlcd_subset[3]-nlcd_subset[2]
    ny=nlcd_subset[1]-nlcd_subset[0]
    x,y=np.meshgrid(x,y)
    x=x.reshape((nx*ny))
    y=y.reshape((nx*ny))

    points=np.concatenate([x[:,np.newaxis],y[:,np.newaxis]],axis=1)

    print("Converting Points from map projection (albers) to lat/lon")
    print(str(data.proj))
    lat,lon=mygis.proj2ll(points=points,proj=data.proj)
    lat=lat.reshape((ny,nx))
    lon=lon.reshape((ny,nx))
    return Bunch(lat=lat,lon=lon)

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

def write(filename,data,geo):
    """write DEM data to a file after subsetting to 'good' data"""
    # goodpoints=np.where(data>0)
    # ymin=np.min(goodpoints[0])
    # ymax=np.max(goodpoints[0])+1
    # xmin=np.min(goodpoints[1])
    # xmax=np.max(goodpoints[1])+1
    xmin=0;xmax=None
    ymin=0;ymax=None
    
    output_data=data[ymin:ymax,xmin:xmax]
    output_lat=geo.lat[ymin:ymax]
    output_lon=geo.lon[xmin:xmax]
    
    lat_info.data=output_lat
    lon_info.data=output_lon
    
    # data_info.subset="[{},{},{},{}]".format(ymin,ymax,xmin,xmax)
    extra_vars=[lat_info,lon_info]
    mygis.write(filename,output_data,varname=data_info.name,dtype=data_info.dtype,dims=data_info.dims,
               attributes=data_info.attributes,extravars=extra_vars,history="wsc/nlcd.py")
    
def load_forest(subset=0):
    """docstring for load_data"""
    global nlcd_subset
    if type(subset)==list:
        nlcd_subset=subset
    elif subset==None:
        pass
    else:
        nlcd_subset=subset_list[subset]
    
    lc_data=mygis.read_tiff(nlcd_file)
    print("setting up geo data")
    lc_geo=geo2lat_lon(lc_data,geo_geographic=True,reverse_lat=True)
    # print("first subset")
    subdata=lc_data.data[nlcd_subset[1]:nlcd_subset[0]:-1,nlcd_subset[2]:nlcd_subset[3]]
    
    lc_forest=np.zeros(subdata.shape)
    lc_forest[np.where((subdata>40)&(subdata<44))]=NLCD_FOREST_TYPE
    
    # these are huge datasets, so help the garbage collector out?
    del subdata, lc_data.data, lc_data
    return Bunch(data=lc_forest, lat=lc_geo.lat,lon=lc_geo.lon)

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
    lc_forest[np.where((subdata>40)&(subdata<44))]=NLCD_FOREST_TYPE
    
    print("Lowering resolution")
    output_data=lower_resolution(lc_forest,lc_geo,geo_data,lc_data.proj)

    print("Writing data file")
    write("snodas_forest_5.nc",output_data,geo_data)
    
if __name__ == '__main__':
    forest_cover()