#!/usr/bin/env python
# (c) Ethan Gutmann gutmann@ucar.edu
from __future__ import print_function
import numpy as np
import mygis

from bunch import Bunch

from wsc import modscag, snodas, wrf, nlcd, regrid
from stat_down import regrid_hi2low as agg

def load_data(load_wrf=True,load_snodas=False):
    """loads either WRF or SNODAS data"""
    if load_snodas:
        return snodas.load("/d2/gutmann/wsc/snodas/SWE_Daily0600UTC_WesternUS_2008.dat",extract_day=121,start_year=2008)
    if load_wrf:
        return wrf.load_base_data(swefile="/d2/gutmann/wsc/wrf/wrfout_d01_2008-05-01_00:00:00",res=2)

def write_data(data):
    """write data to a file, for now just spits out a dummy netcdf file"""
    mygis.write("test_file.nc",data)

def match_grid_nn(inputgrid,outputgrid):
    """nearest neighbor match of an input grid to an output grid"""
    o=outputgrid
    
    lathi,latlo=o.lat,inputgrid.lat
    lonhi,lonlo=o.lon,inputgrid.lon
    winsize=3
    ny,nx=o.data.shape
    ny_l,nx_l=inputgrid.data.shape
    outputdata=np.zeros(o.data.shape)
    
    step=ny/100
    nextoutput=0+step
    for i in range(ny):
        if i>nextoutput:
            print("{0:5.1f}% Completed".format(float(i)/ny*100.0),end="\r")
            nextoutput+=step
        
        dists=(latlo-lathi[i,0])**2+(lonlo-lonhi[i,0])**2
        y,x=np.unravel_index(dists.argmin(),dists.shape)
        for j in range(nx):
            xmax=min(nx_l,x+winsize)
            xmin=max(0,x-winsize)
            ymax=min(ny_l,y+winsize)
            ymin=max(0,y-winsize)
            
            windists=((latlo[ymin:ymax,xmin:xmax]-lathi[i,j])**2+
                      (lonlo[ymin:ymax,xmin:xmax]-lonhi[i,j])**2)
            yoff,xoff=np.unravel_index(windists.argmin(),windists.shape)
            x=xoff+xmin
            y=yoff+ymin
            
            outputdata[i,j]=inputgrid.data[y,x]
    
    return Bunch(data=outputdata,lat=o.lat.copy(),lon=o.lon.copy())
            
def load_hi_res_dem():
    """docstring for load_hi_res_dem"""
    dem_file="/d2/gutmann/wsc/dem/DEM_CO_NHDPlus_1_arc_seconds.nc"
    x0=2500;x1=4000;y0=12000;y1=10500
    x0=0;x1=None;y0=-1;y1=0
    # dem=mygis.read_nc(dem_file,"elev_m").data[y0:y1:-1,x0:x1]
    dem=mygis.read_nc(dem_file,"elev_m").data[::-1,:]
    lat=mygis.read_nc(dem_file,"lat").data[::-1]
    lon=mygis.read_nc(dem_file,"lon").data
    lon,lat=np.meshgrid(lon,lat)
    return Bunch(data=dem,lat=lat,lon=lon)

def write_lc_dem(filename,d1,d2,lat,lon):
    """docstring for write_lc_dem"""
    outputdata=np.zeros((2,d1.shape[0],d1.shape[1]))
    outputdata[0,...]=d1
    outputdata[1,...]=d2
    geo_data=[Bunch(data=lat,dims=("lat","lon"),name="lat"),
              Bunch(data=lon,dims=("lat","lon"),name="lon")]
    
    mygis.write(output_file,outputdata,varname="dem",dims=('lc','lat','lon'),extravars=geo_data)

def create_lc_dem(output_file):
    """docstring for create_lc_dem"""
    dem=load_hi_res_dem()
    lc=nlcd.load_forest(subset=1)
    lc=match_grid_nn(inputgrid=lc,outputgrid=dem)
    
    modscag_grid=modscag.load("/d2/gutmann/wsc/MODSCAG/fsca2008.dat")
    lonlo,latlo=np.meshgrid(modscag_grid.lon,modscag_grid.lat)
    lathi,lonhi = lc.lat,lc.lon
    agg_lut=regrid.agg_lut(lathi,lonhi,latlo,lonlo)
    
    forest_dem_input=np.ma.array(dem.data,mask=(lc.data!=nlcd.NLCD_FOREST_TYPE))
    forest_dem=agg.regrid_hi2low(forest_dem_input,geoLUT=agg_lut,agg_func=np.mean)

    exposed_dem_input=np.ma.array(dem.data,mask=(lc.data!=nlcd.NLCD_EXPOSED_TYPE))
    exposed_dem=agg.regrid_hi2low(exposed_dem_input,geoLUT=agg_lut,agg_func=np.mean)
    
    write_lc_dem(output_file,forest_dem,exposed_dem,latlo,lonlo)
    
    return Bunch(lat=dem.lat,lon=dem.lon,forest=forest_dem,exposed=exposed_dem)
    
    
def load_lc_dem(filename):
    """docstring for load_lc_dem"""
    lat=mygis.read_nc(filename,"lat")
    lon=mygis.read_nc(filename,"lon")
    dem=mygis.read_nc(filename,"dem")
    
    return Bunch(lat=lat,lon=lon,forest=dem[0,...],exposed=dem[1,...])

def load_info():
    """docstring for load_info"""
    lc_dem_file="LC_dem_data"
    try:
        return load_lc_dem(lc_dem_file)
    except:
        return create_lc_dem(lc_dem_file)

def main():
    """convert/downscale snodas and WRF data to modscag grid"""
    modscag_info=load_info()
    data=load_data()
    
    # output_data=rescale(data,modscag_info)
    # write_data(output_data)
    

if __name__ == '__main__':
    main()