#!/usr/bin/env python
import numpy as np
from scipy.signal import medfilt2d
import matplotlib.pyplot as plt

import mygis as io
from bunch import Bunch

import sys
import flushprint
if flushprint.in_ipython():
    pass
else:
    sys.stdout=flushprint.Flushfile(sys.stdout)


def filter_dem(data):
    """Filter a LIDAR DEM to remove some spurious jumps"""
    return medfilt2d(data,5)
    
def calc_veg_height(veg_top,dem):
    veg1=veg_top-dem
    veg1[veg1<0]=0
    veg1[dem>3500]=0
    return medfilt2d(veg1,5)

def calc_snow_depth(snow_data,dem):
    return medfilt2d(snow_data,3)-dem


def decimate(data,factor=10, dfunc=np.mean):
    """Decimate data by a given factor in both dimension
    Uses averaging to decimate by default, could use median, min, max, ...
    """
    sz=data.shape
    
    newsz=(np.floor(data.shape[0]/factor),np.floor(data.shape[1]/factor))
    # the fast method if dfunc support axis= (e.g. most numpy stats: mean,min,max,median,...)
    try:
        data=data[:newsz[0]*factor,:newsz[1]*factor]
        data=data.reshape((newsz[0],factor,newsz[1],factor))
        outputdata=dfunc(dfunc(data,axis=3),axis=1)
    except TypeError:
    # slower but more general way incase dfunc doesn't support axis=  (should throw a TypeError)
        # newsz=(np.round(data.shape[0]/factor),np.round(data.shape[1]/factor))
        outputdata=np.zeros(newsz)
        for i in range(newsz[0]):
            for j in range(newsz[1]):
                # assumes newsz*factor <= len(data)
                outputdata[i,j]=dfunc(data[i*factor:(i+1)*factor,j*factor:(j+1)*factor])
                # assumes len(data) < newsz*factor use this if newsz is computed with round instead of floor
                # if (i<(newsz[0]-1)):
                #     if (j<(newsz[1]-1)):
                #         outputdata[i,j]=dfunc(data[i*factor:(i+1)*factor,j*factor:(j+1)*factor])
                #     else:
                #         outputdata[i,j]=dfunc(data[i*factor:(i+1)*factor,j*factor:])
                # else:
                #     if (j<(newsz[1]-1)):
                #         outputdata[i,j]=dfunc(data[i*factor:,j*factor:(j+1)*factor])
                #     else:
                #         outputdata[i,j]=dfunc(data[i*factor:,j*factor:])
    
    return outputdata
                

def find_veg(vegz):
    """
    Find all locations within a grid of vegetation heights (vegz) that are primarily vegetation
    Decimates by a factor of 10x10 to get a 10m grid
    """
    return vegz>2    

def load(snowoff="compiled_lidar.nc",snowon="snow-on-dem.nc"):
    """Load data from two netCDF files, 
    convert data to snow depth, elevation, and veg height"""
    bare_data=io.read_nc(snowoff,"data").data
    dem=filter_dem(bare_data[0,...])
    veg=calc_veg_height(bare_data[1,...],dem)
    snow_data=io.read_nc(snowon,"Band1").data
    snow=calc_snow_depth(snow_data,dem)
    
    return Bunch(snow=snow,dem=dem,veg=veg,raw=[bare_data,snow_data])

def load_fast(dem="demf",snow="snowdepth",veg="vegmask"):
    """docstring for load_fast"""
    demd=io.read_nc(dem).data
    snowd=io.read_nc(snow).data
    vegd=io.read_nc(veg).data
    return Bunch(veg=vegd,snow=snowd,dem=demd)

def bin_by_elevation(data,dem,mask,dz=50):
    """docstring for bin_by_elevation"""
    minz=dem[dem>100].min()
    maxz=dem[dem<5000].max()
    
    n=np.round((maxz-minz)/dz)
    veg=np.zeros(n)
    vegmed=np.zeros(n)
    vegmin=np.zeros(n)
    vegmax=np.zeros(n)
    exposed=np.zeros(n)
    exposedmed=np.zeros(n)
    exposedmin=np.zeros(n)
    exposedmax=np.zeros(n)
    z=np.arange(minz,maxz+dz,dz)
    
    for i in np.arange(n):
        curexp=np.where((dem>z[i])&(dem<=z[i+1])&(mask==0)&np.isfinite(data)&(data>0)&(data<20))
        curn=len(curexp[0])
        if curn>0:
            exposed[i]=np.mean(data[curexp])
            sorted_data=np.sort(data[curexp])
            exposedmin[i]=sorted_data[int(curn*0.1)]
            exposedmed[i]=sorted_data[int(curn*0.5)]
            exposedmax[i]=sorted_data[int(curn*0.9)]

        curveg=np.where((dem>z[i])&(dem<=z[i+1])&(mask==1)&np.isfinite(data)&(data>0)&(data<20))
        curn=len(curveg[0])
        if curn>0:
            veg[i]=np.mean(data[curveg])
            sorted_data=np.sort(data[curveg])
            vegmin[i]=sorted_data[int(curn*0.1)]
            vegmed[i]=sorted_data[int(curn*0.5)]
            vegmax[i]=sorted_data[int(curn*0.9)]

    veg=np.ma.array(veg,mask=(veg==0))
    vegmin=np.ma.array(vegmin,mask=(vegmin==0))
    vegmed=np.ma.array(vegmed,mask=(vegmed==0))
    vegmax=np.ma.array(vegmax,mask=(vegmax==0))

    return Bunch(z=z[:n],veg=veg,vegmed=vegmed,vegmin=vegmin,vegmax=vegmax,
                exposed=exposed,exposedmed=exposedmed,exposedmin=exposedmin,exposedmax=exposedmax)
    

def main():
    print("Loading data...")
    data=load_fast()
    
    decimation_factor=924 # 924 # for SNODAS or other 30arcsec grid
    print("Decimating Snow")
    snow=decimate(data.snow,decimation_factor)
    print("Decimating Veg")
    veg=decimate(data.veg,decimation_factor,dfunc=np.median)
    # veg_locations=find_veg(veg)
    print("Decimating DEM")
    dem=decimate(data.dem,decimation_factor)
    
    print("Binning")
    banded=bin_by_elevation(snow,dem,veg,dz=100)
    
    print("Plotting")
    plt.plot(banded.z,banded.exposed,label="Exposed",color="b",linewidth=2)
    plt.plot(banded.z,banded.exposedmed,"--",label="Exp. Median",color="b",linewidth=2)
    plt.bar([0],[1],color="skyblue",edgecolor="black",label="Exp. 10-90%")
    plt.plot(banded.z,banded.veg,label="Vegetation",color="g",linewidth=2)
    plt.plot(banded.z,banded.vegmed,"--",label="Veg. Median",color="g",linewidth=2)
    plt.bar([0],[1],color="lightgreen",edgecolor="black",label="Veg. 10-90%")
    
    plt.plot(banded.z,banded.vegmin,color="black")
    plt.plot(banded.z,banded.vegmax,color="black")
    plt.plot(banded.z,banded.exposedmin,color="black")
    plt.plot(banded.z,banded.exposedmax,color="black")

    plt.fill_between(banded.z,banded.exposedmin,banded.exposedmax,
                        color="skyblue",edgecolor="black")
    plt.fill_between(banded.z,banded.vegmin,banded.vegmax,
                        color="lightgreen",alpha=0.5,edgecolor="black")
                        

    plt.legend(loc=2)
    plt.xlim(2500,3800)
    plt.ylim(0,6)
    plt.ylabel("Snow Depth [m]")
    plt.xlabel("Elevation [m]")
    plt.title("Lidar snow depth at Niwot")
    plt.savefig("lidar_by_elev.png")
    
    
    
    

if __name__ == '__main__':
    main()