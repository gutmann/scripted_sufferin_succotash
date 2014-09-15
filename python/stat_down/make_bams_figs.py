#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

from stat_down import map_vis
import mygis
from bunch import Bunch

def summarize_precip(current_data, future_data):

    cur=current_data[0][-1,:,:] - current_data[0][0,:,:]
    fut=future_data[0][-1,:,:]  - future_data[0][0,:,:]
    for i in range(1,len(current_data)):
        cur+=(current_data[i][-1,:,:] - current_data[i][0,:,:])
        fut+=(future_data[i][-1,:,:]  - future_data[i][0,:,:])
    
    cur/=len(current_data)
    fut/=len(current_data)
    
    return Bunch(current=cur,future=fut)

def summarize_temp(current_data, future_data):
    """docstring for summarize_snow"""
    cur=current_data[0].mean(axis=0)
    fut=future_data[0].mean(axis=0)
    for i in range(1,len(current_data)):
        cur+=current_data[i].mean(axis=0)
        fut+=future_data[i].mean(axis=0)
    cur/=len(current_data)
    fut/=len(current_data)
    
    return Bunch(current=cur,future=fut)


def summarize_snow(current_data, future_data):
    """docstring for summarize_snow"""
    cur=current_data[0][-1,:,:]
    fut=future_data[0][-1,:,:]
    for i in range(1,len(current_data)):
        cur+=current_data[i][-1,:,:]
        fut+=future_data[i][-1,:,:]
    cur/=len(current_data)
    fut/=len(current_data)
    
    return Bunch(current=cur,future=fut)


def load_wrf_march_data(varname):
    """docstring for load_wrf_march_data"""
    # current_data=mygis.read_files("NARR*_"+varname+"_*2004*.nc",varname)
    # future_data=mygis.read_files("PGW*_"+varname+"_*2004*.nc",varname)
    current_data=mygis.read_files("NARR*_"+varname+"_*.nc",varname)
    future_data=mygis.read_files("PGW*_"+varname+"_*.nc",varname)
    
    if varname=="SNOW":
        return summarize_snow(current_data,future_data)
    
    if varname=="RAINNC":
        return summarize_precip(current_data,future_data)
    
    if varname=="T2":
        return summarize_temp(current_data,future_data)
        
def load_wrf_elevation():
    """docstring for load_wrf_elevation"""
    wrf_baseline_file="/glade/u/home/gutmann/scratch/icar/baseline_info/4km_headwaters_subset.nc"
    return mygis.read_nc(wrf_baseline_file,"HGT").data[0,...]
        

def make_plots(data,elevation):
    """docstring for make_plots"""
    print("plotting data")
    plt.figure(figsize=(12,8))
    plt.subplot(2,2,1)
    vname="T2"
    map_vis.vis(data[vname].future - data[vname].current,geo="subset",proj="lcc",latstep=2.0,lonstep=5.0,
        cmap=plt.cm.jet,title="Temperature",cbar_label="[degrees C]")
        # cmap=plt.cm.seismic,clim=(-5,5),title="Temperature",cbar_label="[degrees C]")

    plt.subplot(2,2,2)
    vname="RAINNC"
    map_vis.vis(data[vname].future - data[vname].current,geo="subset",proj="lcc",latstep=2.0,lonstep=5.0,
        cmap=plt.cm.seismic_r,clim=(-50,50),title="Precipitation",cbar_label="[millimeters]")
    
    plt.subplot(2,2,3)
    map_vis.vis(elevation,geo="subset",proj="lcc",latstep=2.0,lonstep=5.0,
        cmap=plt.cm.terrain,clim=(-400,3800),title="Terrain",cbar_label="[meters]")
    
    plt.subplot(2,2,4)
    vname="SNOW"
    map_vis.vis(data[vname].future - data[vname].current,geo="subset",proj="lcc",latstep=2.0,lonstep=5.0,
        cmap=plt.cm.seismic_r,clim=(-100,100),title="SWE",cbar_label="[millimeters]")
    
    plt.savefig("wrf_march_fig.png",dpi=200)
    

def wrf_march_fig():
    """docstring for wrf_march_fig"""
    variables=["T2","RAINNC","SNOW"]
    data=dict()
    print("loading elevation")
    elevation=load_wrf_elevation()
    
    for v in variables:
        print("loading "+v)
        data[v]=load_wrf_march_data(v)
    
    make_plots(data,elevation)

def main():
    """docstring for main"""
    wrf_march_fig()

if __name__ == '__main__':
    main()