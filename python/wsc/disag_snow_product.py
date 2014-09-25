#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from wsc import snodas,wrf,modscag

from bunch import Bunch 
import mygis

filenames=dict(wrf="wrf_high_res.nc",snodas="snodas_high_res.nc")

def load_wrf():
    """docstring for load_wrf"""
    wrffilename="wrf/SWE_daily.nc"
    wrfgeofile="wrf/4km_wrf_output.nc"
    
    data=mygis.read_nc(wrffilename,"SNOW").data[-365:,...].max(axis=0)/1000.0
    lat=mygis.read_nc(wrfgeofile,"XLAT").data[0,...]
    lon=mygis.read_nc(wrfgeofile,"XLONG").data[0,...]
    hgt=mygis.read_nc(wrfgeofile,"HGT").data[0,...]
    
    return Bunch(snow=data,lat=lat,lon=lon,z=hgt)
    
def load_snodas():
    """docstring for load_snodas"""
    snodasfilename="snodas/SWE_Daily0600UTC_WesternUS_2008.dat"
    data=snodas.load(snodasfilename,startyear=2008,fill=True)
    if len(data.lon.shape)<2:
        lon,lat=np.meshgrid(data.lon,data.lat)
    else:
        lon,lat=(data.lon,data.lat)
    return Bunch(snow=data.data.max(axis=0),lat=lat,lon=lon,z=data.dem)

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
        print("{:6.2f}%".format(100.0*i/float(ny)),end="\r")
        for j in range(nx):
            cury,curx=find_point(grid.lat[i,j],grid.lon[i,j],data.lat,data.lon,lasty,lastx)
            miny=max(cury-window_size,0)
            maxy=min(cury+window_size+1,ny_low)
            minx=max(curx-window_size,0)
            maxx=min(curx+window_size+1,nx_low)
            slope,offset=np.polyfit(data.z[miny:maxy,minx:maxx].flat,data.snow[miny:maxy,minx:maxx].flat,deg=1)
            exposed[i,j]=grid.exposed.z[i,j]*slope+offset
            forest[i,j] =grid.forest.z[i,j] *slope+offset
            
            lastx,lasty=curx,cury
    print("")
    return Bunch(exposed=exposed,forest=forest)

def write_output(data,filename):
    """docstring for write_output"""
    mygis.write("forest_"+filename,data.forest)
    mygis.write("exposed_"+filename,data.exposed)
    

def main(method="snodas"):
    """docstring for main"""
    print("Loading low-res data")
    if method=="wrf":
        data=load_wrf()
    elif method=="snodas":
        data=load_snodas()
    
    print("Loading high-res data")
    grid=load_high_res_grid()
    
    print("Downscaling")
    high_res=downscale_snow(data,grid)
    
    print("Writing output")
    write_output(high_res,filenames[method])

if __name__ == '__main__':
    main()