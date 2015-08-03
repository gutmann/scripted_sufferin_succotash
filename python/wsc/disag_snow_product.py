#!/usr/bin/env python
"""
usage: disag_snow_product.py [-h] [-m METHOD] [-y YEAR] [-v] [--verbose]

Disaggregate WRF or SNODAS SWE data to the MODSCAG grid

optional arguments:
  -h, --help            show this help message and exit
  -m METHOD, --method METHOD
                        dataset to process [wrf,snodas] (default: wrf)
  -y YEAR, --year YEAR  Year to process [2004-2008] (default: 2004)
  -v, --version         show program's version number and exit
  --verbose             verbose output (default: False)
"""
from __future__ import print_function

import sys
import os
import traceback
import argparse

import numpy as np
from wsc import snodas,wrf,modscag

from bunch import Bunch 
import mygis

filenames=dict(wrf="wrf_high_res.nc",snodas="snodas_high_res.nc")

def load_wrf(year=2004, load_may=False):
    """docstring for load_wrf"""
    wrffilename="wrf/SWE_daily.nc"
    wrfgeofile="wrf/4km_wrf_output.nc"
    startpt=int((year-2001)*365.25)
    data=mygis.read_nc(wrffilename,"SNOW").data[startpt:startpt+365,...]
    if load_may:
        data=data[212]/1000.0
    else:
        data=data.max(axis=0)/1000.0
    lat=mygis.read_nc(wrfgeofile,"XLAT").data[0,...]
    lon=mygis.read_nc(wrfgeofile,"XLONG").data[0,...]
    hgt=mygis.read_nc(wrfgeofile,"HGT").data[0,...]
    
    return Bunch(snow=data,lat=lat,lon=lon,z=hgt)
    
def load_snodas(year=2004, load_may=False):
    """docstring for load_snodas"""
    snodasfilename="snodas/SWE_Daily0600UTC_WesternUS_{year}.dat".format(year=year)
    data=snodas.load(snodasfilename,startyear=year,fill=True)
    if len(data.lon.shape)<2:
        lon,lat=np.meshgrid(data.lon,data.lat)
    else:
        lon,lat=(data.lon,data.lat)
    
    if load_may:
        snowdata=data.data[120]
    else:
        snowdata=data.data.max(axis=0)
    
    return Bunch(snow=snowdata,lat=lat,lon=lon,z=data.dem)

def load_high_res_grid():
    """docstring for load_high_res_grid"""
    #Load data
    modscag_file="MODSCAG/fsca2008.dat"
    exp_dem=mygis.read_nc("exposed_dem.nc").data
    for_dem=mygis.read_nc("forest_dem.nc").data
    exp_n=mygis.read_nc("exposed_n.nc").data
    for_n=mygis.read_nc("forest_n.nc").data
    grid=modscag.load(modscag_file)
    
    # find valid subset
    xn=(exp_n+for_n).sum(axis=0)
    yn=(exp_n+for_n).sum(axis=1)
    goodx=np.where(xn>0)
    goody=np.where(yn>0)
    minx=goodx[0][0]
    maxx=goodx[0][-1]
    miny=goody[0][0]
    maxy=goody[0][-1]
    
    # subset data
    internal_exposed_n=exp_n[miny:maxy,minx:maxx].copy()
    internal_exposed_n[internal_exposed_n==0]=1
    internal_forest_n=for_n[miny:maxy,minx:maxx].copy()
    internal_forest_n[internal_forest_n==0]=1
    
    lon,lat=np.meshgrid(grid.lon[minx:maxx],grid.lat[miny:maxy])
    
    
    return Bunch(exposed=Bunch(z=exp_dem[miny:maxy,minx:maxx]/internal_exposed_n, n=exp_n[miny:maxy,minx:maxx]),
                 forest =Bunch(z=for_dem[miny:maxy,minx:maxx]/internal_forest_n,  n=for_n[miny:maxy,minx:maxx]),
                 lat=lat,lon=lon)

def find_point(lat_pt,lon_pt,lat,lon,lasty,lastx):
    """docstring for find_point"""
    
    if lasty<0:
        dists=(lat-lat_pt)**2 + (lon-lon_pt)**2
        return np.unravel_index(np.argmin(dists),dists.shape)
    
    else:
        window_size=2
        ny,nx=lat.shape
        
        miny=max(lasty-window_size,0)
        maxy=min(lasty+window_size+1,ny)
        minx=max(lastx-window_size,0)
        maxx=min(lastx+window_size+1,nx)
        dists=(lat[miny:maxy,minx:maxx]-lat_pt)**2 + (lon[miny:maxy,minx:maxx]-lon_pt)**2
        cury,curx=np.unravel_index(np.argmin(dists),dists.shape)
        
        return cury+miny,curx+minx
        
    
def downscale_snow(data,grid):
    """docstring for downscale_snow"""
    ny,nx=grid.lat.shape
    ny_low,nx_low=data.lat.shape
    
    lasty,lastx=-1,-1
    window_size=1
    exposed=np.zeros((ny,nx))
    forest=np.zeros((ny,nx))
    
    for i in range(ny):
        lasty,lastx=-1,-1
        if verbose: print("{:6.2f}%".format(100.0*i/float(ny)),end="\r")
        for j in range(nx):
            cury,curx=find_point(grid.lat[i,j],grid.lon[i,j],data.lat,data.lon,lasty,lastx)
            miny=max(cury-window_size,0)
            maxy=min(cury+window_size+1,ny_low)
            minx=max(curx-window_size,0)
            maxx=min(curx+window_size+1,nx_low)
            
            slope,offset=np.polyfit(data.z[miny:maxy,minx:maxx].flat,data.snow[miny:maxy,minx:maxx].flat,deg=1)
            if (slope<0):
                slope=0
                offset=data.snow[cury,curx]
                
            exposed[i,j]=grid.exposed.z[i,j]*slope+offset
            forest[i,j] =grid.forest.z[i,j] *slope+offset
            
            lastx,lasty=curx,cury
    if verbose: print("")
    return Bunch(exposed=exposed,forest=forest)

def write_output(data,filename,year,may):
    """docstring for write_output"""
    if may:
        filename="May1st_"+filename
    mygis.write("forest_{year}_{file}".format(year=year,file=filename), data.forest)
    mygis.write("exposed_{year}_{file}".format(year=year,file=filename),data.exposed)
    

def main(method="snodas", year=2004, may=False):
    """docstring for main"""
    if verbose: print("Loading low-res data")
    if method=="wrf":
        data=load_wrf(year,load_may=may)
    elif method=="snodas":
        data=load_snodas(year,load_may=may)
    
    if verbose: print("Loading high-res data")
    grid=load_high_res_grid()
    
    if verbose: print("Downscaling")
    high_res=downscale_snow(data,grid)
    
    if verbose: print("Writing output")
    write_output(high_res,filenames[method],year,may)


if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Disaggregate WRF or SNODAS to MODSCAG grid. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-m','--method',dest="method", action='store', default="wrf", help="dataset to process [wrf,snodas]")
        parser.add_argument('-y','--year',  dest="year",   action='store', default=2004,  help="Year to process [2004-2008]", type=int)
        parser.add_argument ('--may', action='store_true',
                default=False, help='Load May1st SWE instead of seasonal maximum', dest='may')
        parser.add_argument('-v', '--version',action='version',
                version='disag_snow_product 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
        
        global verbose
        verbose=args.verbose

        exit_code = main(method=args.method,year=args.year, may=args.may)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
