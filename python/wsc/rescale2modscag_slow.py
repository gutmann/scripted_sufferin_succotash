#!/usr/bin/env python
# (c) Ethan Gutmann gutmann@ucar.edu
from __future__ import print_function
import sys
import numpy as np
import mygis

from bunch import Bunch

from wsc import modscag, snodas, wrf, nlcd

dem_file="/d2/gutmann/wsc/dem/DEM_CO_NHDPlus_1_arc_seconds.nc"
modscag_file="/d2/gutmann/wsc/MODSCAG/fsca2008.dat"
nlcd_file="/d2/gutmann/nldc/NLCD2006_landcover.tif"
nlcd_subset=[25000,75000,30000,63000]
DEBUG=False

def load_hi_res_dem():
    """docstring for load_hi_res_dem"""
    demf=mygis.read_nc(dem_file,"elev_m",returnNCvar=True)
    latf=mygis.read_nc(dem_file,"lat",returnNCvar=True)
    lonf=mygis.read_nc(dem_file,"lon",returnNCvar=True)
    if DEBUG:
        x0=2500;x1=4000;y0=12000;y1=10500
    # x0=0;x1=None;y0=-1;y1=0
        dem=demf.data[y0:y1:-1,x0:x1]
        lat=latf.data[y0:y1:-1]
        lon=lonf.data[x0:x1]
    else:
        dem=demf.data[::-1,:]
        lat=latf.data[::-1]
        lon=lonf.data[:]
    demf.ncfile.close()
    latf.ncfile.close()
    lonf.ncfile.close()
    
    dx=lon[1]-lon[0]
    dy=lat[1]-lat[0]
    nx=len(lon)
    ny=len(lat)
    lon,lat=np.meshgrid(lon,lat)
    
    return Bunch(data=dem,lat=lat,lon=lon,
                 startx=lon[0,0]-dx/2,starty=lat[0,0]-dy/2,
                 dx=dx,dy=dy,nx=nx,ny=ny)

def load_forest():
    """docstring for load_forest"""
    
    lc_data=mygis.read_tiff(nlcd_file)
    subdata=lc_data.data[nlcd_subset[1]:nlcd_subset[0]:-1,nlcd_subset[2]:nlcd_subset[3]]
    left=lc_data.geo[0]
    top =lc_data.geo[3]
    dx  =lc_data.geo[1]
    dy  =lc_data.geo[5]
    
    nx=lc_data.data.shape[1]
    ny=lc_data.data.shape[0]
    
    x=np.arange(left,left+dx*nx,dx)[nlcd_subset[2]:nlcd_subset[3]].astype(np.int32)
    y=np.arange(top, top +dy*ny,dy)[nlcd_subset[1]:nlcd_subset[0]:-1].astype(np.int32)
    nx=len(x)
    ny=len(y)
    dy= -1*dy
    lc_forest=np.zeros(subdata.shape)
    if (nlcd.NLCD_EXPOSED_TYPE!=0):
        lc_forest[:]=nlcd.NLCD_EXPOSED_TYPE
    lc_forest[np.where((subdata>40)&(subdata<44))]=nlcd.NLCD_FOREST_TYPE
    
    proj=lc_data.proj
    # these are huge datasets, so help the garbage collector out?
    del subdata, lc_data.data, lc_data
    return Bunch(data=lc_forest, x=x,y=y,proj=proj, nx=nx,ny=ny,startx=x[0],starty=y[0],dx=dx,dy=dy)

def write_output(data):
    """docstring for write_output"""
    for k in data.keys():
        mygis.write(k+".nc",data[k])

def find_point(grid,findx,findy):
    """docstring for get_point"""
    outx=(findx-grid.startx)/grid.dx
    outy=(findy-grid.starty)/grid.dy
    
    if ((outx>=0) and (outx<grid.nx) and (outy>=0) and (outy<grid.ny)):
        return outx,outy
    else:
        return None,None
            

def main():
    """docstring for main"""
    print("loading dem")
    dem=load_hi_res_dem()
    print("loading LC")
    lc=load_forest()
    print(lc.startx,lc.starty,lc.dx,lc.dy)
    print(lc.proj)
    
    print("loading modscag")
    modscag_grid=modscag.load(modscag_file)
    nx=len(modscag_grid.lon)
    ny=len(modscag_grid.lat)
    ms_gridinfo=Bunch(startx=modscag_grid.lon[0],starty=modscag_grid.lat[0],
                 nx=nx,ny=ny,dx=modscag_grid.lon[1]-modscag_grid.lon[0],
                 dy=modscag_grid.lat[1]-modscag_grid.lat[0])
    print(ms_gridinfo)
    outputdata=Bunch(forest_dem=np.zeros((ny,nx)),exposed_dem=np.zeros((ny,nx)),
                     forest_n=np.zeros((ny,nx)),  exposed_n=np.zeros((ny,nx)))
    
    nx=dem.data.shape[1]
    ny=dem.data.shape[0]
        
    for i in range(ny):
        print("{0:5} of {1:5} = {2:5.1f}%".format(i,ny,100*i/float(ny)),end="\r")
        sys.stdout.flush()
        for j in range(nx):
            curlat=dem.lat[i,j]
            curlon=dem.lon[i,j]
            msx,msy=find_point(ms_gridinfo,curlon,curlat)
            
            if msx:
                # print("found:",msx,msy,curlon,curlat,modscag_grid.lon[msx],modscag_grid.lat[msy])
                curx,cury=mygis.ll2proj(lon=curlon,lat=curlat,proj=lc.proj)
                # print(curx,cury,curx-lc.startx,lc.dx,cury-lc.starty,lc.dy)
                lcx,lcy=find_point(lc,curx,cury)
                if lcx:
                    # print("LC found:",lcx,lcy,lc.x[lcx],lc.y[lcy],curx,cury)
                    if (lc.data[lcy,lcx]==nlcd.NLCD_FOREST_TYPE):
                        outputdata.forest_dem[msy,msx]+=dem.data[i,j]
                        outputdata.forest_n[msy,msx]+=1
                    else:
                        outputdata.exposed_dem[msy,msx]+=dem.data[i,j]
                        outputdata.exposed_n[msy,msx]+=1
    print("")
    print("Writing data")
    write_output(outputdata)
    
    
            


if __name__ == '__main__':
    main()