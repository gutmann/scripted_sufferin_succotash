#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys

import mygis
from bunch import Bunch

import wsc.lidar as lidar
import wsc.wrf as wrf
import wsc.snodas as snodas

def zband(data,valid,weights=None):
    curn=len(valid[0])
    if curn>0:
        if weights!=None:
            output_vol=np.sum(data[valid]*weights[valid])
            output=output_vol/np.sum(weights[valid])
        else:
            output=np.mean(data[valid])
            output_vol=output*(curn)
        if curn>2:
            sorted_data=np.sort(data[valid])
            output_min=sorted_data[int(curn*0.1)]
            output_med=sorted_data[int(np.floor(curn*0.5))]
            output_max=sorted_data[int(np.floor(curn*0.9))]
        else:
            output_min=np.nan
            output_med=np.nan
            output_max=np.nan
    else:
        output=np.nan
        output_vol=np.nan
        output_min=np.nan
        output_med=np.nan
        output_max=np.nan
    
    return Bunch(mean=output,min=output_min,median=output_med,max=output_max,volume=output_vol)

def bin_by_elevation(data,dem,mask,dz=100,weights=None,dx=1):
    """docstring for bin_by_elevation"""
    minz=dem[dem>100].min()
    maxz=dem[dem<5000].max()
    
    n=np.round((maxz-minz)/dz)
    z=np.arange(minz,maxz+dz,dz)
    
    veg=np.zeros(n)
    vegvol=np.zeros(n)
    vegmed=np.zeros(n)
    vegmin=np.zeros(n)
    vegmax=np.zeros(n)
    vegarea=np.zeros(n)
    exposed=np.zeros(n)
    exposedvol=np.zeros(n)
    exposedmed=np.zeros(n)
    exposedmin=np.zeros(n)
    exposedmax=np.zeros(n)
    exposedarea=np.zeros(n)

    
    for i in np.arange(n):
        curexp=np.where((dem>z[i])&(dem<=z[i+1])&(mask==2)&np.isfinite(data)&(data>0)&(data<20))
        band=zband(data,curexp,weights=weights)
        exposed[i]=band.mean
        exposedvol[i]=band.volume*dx**2/dz
        exposedmin[i]=band.min
        exposedmed[i]=band.median
        exposedmax[i]=band.max
        exposedarea[i]=len(curexp[0])*dx**2 / dz
        
        curveg=np.where((dem>z[i])&(dem<=z[i+1])&(mask==1)&np.isfinite(data)&(data>0)&(data<20))
        band=zband(data,curveg,weights=weights)
        veg[i]=band.mean
        vegvol[i]=band.volume*dx**2/dz
        vegmin[i]=band.min
        vegmed[i]=band.median
        vegmax[i]=band.max
        vegarea[i]=len(curveg[0])*dx**2 / dz

    exposed   =np.ma.array(exposed,   mask=(~np.isfinite(exposed)))
    exposedvol=np.ma.array(exposedvol,mask=(~np.isfinite(exposedvol)))
    exposedmin=np.ma.array(exposedmin,mask=(~np.isfinite(exposedmin)))
    exposedmed=np.ma.array(exposedmed,mask=(~np.isfinite(exposedmed)))
    exposedmax=np.ma.array(exposedmax,mask=(~np.isfinite(exposedmax)))
    
    veg   =np.ma.array(veg,   mask=(~np.isfinite(veg)))
    vegvol=np.ma.array(vegvol,mask=(~np.isfinite(vegvol)))
    vegmin=np.ma.array(vegmin,mask=(~np.isfinite(vegmin)))
    vegmed=np.ma.array(vegmed,mask=(~np.isfinite(vegmed)))
    vegmax=np.ma.array(vegmax,mask=(~np.isfinite(vegmax)))

    return Bunch(z=z[:n]+(dz/2),veg=veg,  vegvol=vegvol, vegmed=vegmed,         vegmin=vegmin,         vegmax=vegmax,        vegarea=vegarea,
                exposed=exposed,  exposedvol=exposedvol, exposedmed=exposedmed, exposedmin=exposedmin, exposedmax=exposedmax,exposedarea=exposedarea)

def plot_elevation_bands(banded,outputfile="swe_by_z_lc.png",title="SWE by elevation and Land Cover",
                        ylim=(0,1),ylabel="Snow Water Equivalent [m]", 
                        xlim=(2500,4000),xlabel="Elevation [m]"):
    """docstring for plot_elevation_bands"""
    import matplotlib.pyplot as plt

    print("Plotting")
    plt.clf();
    plt.plot(banded.z,banded.exposed,label="Exposed",color="b",linewidth=2)
    plt.plot(banded.z,banded.exposedmed,"--",label="Exp. Median",color="b",linewidth=2)
    plt.bar([0],[1],color="skyblue",alpha=0.5,edgecolor="black",label="Exp. 10-90%") #note, this is just to add the label to the legend
    plt.plot(banded.z,banded.veg,label="Vegetation",color="g",linewidth=2)
    plt.plot(banded.z,banded.vegmed,"--",label="Veg. Median",color="g",linewidth=2)
    plt.bar([0],[1],color="lightgreen",alpha=0.5,edgecolor="black",label="Veg. 10-90%") # add a label to the legend
    
    plt.plot(banded.z,banded.vegmin,color="black")
    plt.plot(banded.z,banded.vegmax,color="black")
    plt.plot(banded.z,banded.exposedmin,color="black")
    plt.plot(banded.z,banded.exposedmax,color="black")

    plt.fill_between(banded.z,banded.exposedmin,banded.exposedmax,
                        color="skyblue",edgecolor="black")
    plt.fill_between(banded.z,banded.vegmin,banded.vegmax,
                        color="lightgreen",alpha=0.5,edgecolor="black")
                        

    plt.legend(loc=2)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    plt.savefig(outputfile)



def gen_weights(hlat,hlon,llat,llon,mask=None,verbose=False):
    """compute weights for a nearestneighbor sampling between low and high resolution grids
    llat,llon specify the lat and longitude grids on the low resolution grid 
    hlat,hlon specify the lat and longitude grids on the high resolution grid 

    returns a low-res grid with values of n high res points per grid cell
    """
    if len(llat.shape)==1:
        llon,llat=np.meshgrid(llon,llat)
    if len(hlat.shape)==1:
        hlon,hlat=np.meshgrid(hlon,hlat)
    
    output=np.zeros(llon.shape)
    if mask==None:
        mask=np.ones(hlat.shape,dtype=bool)

    search=2
    if verbose:
        print("Total={}".format(hlat.shape[0]))
        
    for i in range(hlat.shape[0]):
        dists=(llat-hlat[i,0])**2 + (llon-hlon[i,0])**2
        lastx,lasty=np.unravel_index(np.argmin(dists),dists.shape)
        if verbose:
            print(i,end=" ")
            sys.stdout.flush()
        for j in range(hlat.shape[1]):
            sx,ex=lastx-search,lastx+search
            sy,ey=lasty-search,lasty+search
            dists=(llat[sx:ex,sy:ey]-hlat[i,j])**2 + (llon[sx:ex,sy:ey]-hlon[i,j])**2
            subx,suby=np.unravel_index(np.argmin(dists),dists.shape)
            curx=subx+sx
            cury=suby+sy
            
            if mask[i,j]:
                output[curx,cury]+=1
                    
            lastx=curx
            lasty=cury
    return output
            
def find_bounds(data):
    """docstring for find_bounds"""
    xaxis=data.sum(axis=0)
    yaxis=data.sum(axis=1)
    xpos=np.where(xaxis>0)[0]
    ypos=np.where(yaxis>0)[0]
    
    return [ypos[0],ypos[-1],xpos[0],xpos[-1]]
    
def plot_volumes(wrf,sno,lidar=None):
    """docstring for plot_volumes"""
    import matplotlib.pyplot as plt
    plt.clf()
    # wrf/sno can be lists (each list element is a different year) or single year
    if type(wrf)==list:
        # add labels for the legend
        # WRF has no years overlapping the lidar, so don't plot anything bold
        plt.plot(wrf[0].z,wrf[0].exposedvol*-1,"--",color="red",label="WRF Exposed")
        plt.plot(wrf[0].z,wrf[0].vegvol*-1,"-", color="red",label="WRF Forest")
        for w in wrf:
            plt.plot(w.z,w.exposedvol,"--",color="salmon")
            plt.plot(w.z,w.vegvol,    "-", color="salmon")
    else:
        plt.plot(wrf.z,wrf.exposedvol,"--",color="red",label="WRF Exposed")
        plt.plot(wrf.z,wrf.vegvol,    "-", color="red",label="WRF Forest")
        
    if type(sno)==list:
        # add labels for the legend
        # assume that the corresponding SNODAS year has been placed at the beginning of the list and plot it bolder
        plt.plot(sno[0].z,sno[0].exposedvol,"--",linewidth=2,color="green",label="SNODAS Exposed")
        plt.plot(sno[0].z,sno[0].vegvol,    "-", linewidth=2,color="green",label="SNODAS Forest")
        for s in sno[1:]:
            plt.plot(s.z,s.exposedvol,"--",color="lightgreen")
            plt.plot(s.z,s.vegvol,    "-", color="lightgreen")
    else:
        plt.plot(sno.z,sno.exposedvol,"--",color="green",label="SNODAS Exposed")
        plt.plot(sno.z,sno.vegvol,    "-", color="green",label="SNODAS Forest")

    if lidar!=None:
        plt.plot(lidar.z,lidar.exposedvol,"--",linewidth=2,color="blue",label="lidar Exposed")
        plt.plot(lidar.z,lidar.vegvol,    "-", linewidth=2,color="blue",label="lidar Forest")
    
    plt.ylim(0,200000)
    plt.xlim(2500,4000)
    plt.ylabel("Water Volume per meter (z) [m$^3$/m]")
    plt.xlabel("Elevation [m]")
    plt.legend(loc=2)
    plt.savefig("SWE_volume_by_z.png")


def main():
    """Compare WRF, SNODAS, and Lidar data on the lidar domain"""
    snowdensity=0.35 #from May 1 2010 SNOTEL (2011,2013 were similar, 2014 was 0.4), at the saddle in May 1 2010 it was 0.4
    snodasyears=[2010,2004,2005]
    wdata=[wrf.load("wrf/SWE_daily.nc",extractday=212+5+int(np.round(365.25*year))) for year in [3,4]]
    wdata.extend([wrf.load("wrf/SWE_daily.nc",extractday=212+20+int(np.round(365.25*year))) for year in [3,4]])
    print(len(wdata))
    sdata=[snodas.load("snodas/SWE_Daily0600UTC_WesternUS_{}.dat".format(year),extractday=125) for year in snodasyears]
    sdata.extend([snodas.load("snodas/SWE_Daily0600UTC_WesternUS_{}.dat".format(year),extractday=140) for year in snodasyears])
    print(len(sdata))
    # sdata=[snodas.load("snodas/SWE_Daily0600UTC_WesternUS_{}.dat".format(year),extractday=120) for year in range(2004,2013)]
    # sdata.insert(0,sdata.pop(6)) #move year 2010 to the begining of the list
    ldata=lidar.load_fast(loc="lidar/",geofile="snow-on-dem.nc",decimation_factor=10)
    
    print("Calculating WRF weights")
    try:
        wrfweights=mygis.read_nc("wrf2lidar_weights.nc").data
    except:
        wrfweights   =gen_weights(ldata.lat,ldata.lon,wdata[0].lat,wdata[0].lon,mask=(ldata.dem>1500))
        mygis.write("wrf2lidar_weights.nc",wrfweights)
        
    # wrfbounds    =find_bounds(wrfweights)
    print("Calculating SNODAS weights")
    try:
        snodasweights=mygis.read_nc("snodas2lidar_weights.nc").data
    except:
        snodasweights=gen_weights(ldata.lat,ldata.lon,sdata[0].lat,sdata[0].lon,mask=(ldata.dem>1500))
        mygis.write("snodas2lidar_weights.nc",snodasweights)
        
    # snodasbounds =find_bounds(snodasweights)
    
    wdata[0].lc[wrfweights==0]=0
    sdata[0].lc[snodasweights==0]=0

    print("Binning by elevations...")
    #dx=4000) #note use dx=lidar_dx because weights are lidar gridcells...
    wrfbyz=[bin_by_elevation(w.data,w.dem,wdata[0].lc,weights=wrfweights,dz=200,dx=10) for w in wdata]
    print("Binning by elevations...")
    snodasbyz=[bin_by_elevation(s.data,sdata[0].dem,sdata[0].lc,weights=snodasweights,dz=150,dx=10) for s in sdata]#dx=926)
    print("Binning by elevations...")
    lidarbyz=bin_by_elevation(ldata.data*snowdensity,ldata.dem,ldata.lc,dz=100,dx=10)
    print("Plotting")
    plot_volumes(wrfbyz,snodasbyz,lidarbyz)

    snodasyears=[2010,2004,2005,2010.2,2004.2,2005.2]
    for i in range(len(snodasbyz)):
        plot_elevation_bands(snodasbyz[i],outputfile="SNODAS_swe_by_z_{}.png".format(snodasyears[i]),title="SNODAS SWE {}".format(snodasyears[i]))
    
    
if __name__ == '__main__':
    main()