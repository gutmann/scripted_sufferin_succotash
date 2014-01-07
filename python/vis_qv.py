#!/usr/bin/env python
from __future__ import print_function
import sys
from glob import glob
import datetime
import date_fun

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from stat_down import myio as io

def load_data(filesearch):
    """docstring for load_data"""
    data=io.read_files(filesearch,"qv")
    return data

def map_vis(data,geo=[],title="",vmin=None,vmax=None,cmap=None,showcolorbar=True,
            latlabels=[1,0,0,0],lonlabels=[0,0,0,1],cbar_label=None):
    """Plot a map of data using the bounds in geo=[lllat,urlat,lllon,urlon]
    
    Optionally specify a map title, min and max value and colormap
    """
    # geo=[35,43,-113,-101]
    geo=[25.125,52.875,-124.75,-67]
    nx=data.shape[1]
    ny=data.shape[0]
    dx=12000.
    m = Basemap(width=nx*dx,height=ny*dx,
                rsphere=(6378137.00,6356752.3142),\
                resolution='l',area_thresh=10000.,projection='lcc',\
                lat_1=28.,lat_2=50.,lat_0=39.7,lon_0=-98.0)

    mapimg=m.imshow(data,vmin=vmin,vmax=vmax,cmap=cmap)
    m.drawparallels(np.arange(25,55,5.),labels=latlabels,dashes=[1,4])
    m.drawmeridians(np.arange(-120,-65,10.),labels=lonlabels,dashes=[1,4])
    m.drawstates(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    
    if showcolorbar:
        cb=m.colorbar()
        if cbar_label!=None:
            cb.set_label(cbar_label)
        
    plt.title(title,x=0.05,ha="left")

    
def plt_data(data,files,dates):
    minval=0.0005
    period=365.25
    offset=20
    max_range=0.005
    max_mean=0.0095
    
    for i in range(381,len(files)):
        plt.clf();
        data=io.read_nc(files[i],"qv").data
        doy=date_fun.datetime2mjd([dates[i]])-date_fun.date2mjd(dates[i].year,1,1,0,0)
        print(dates[i])
        maxval=np.cos((doy-offset-period/2)/period*2*np.pi)*max_range+max_mean
        # print(maxval)
        map_vis(data[:,0,:],vmin=minval,vmax=maxval,cmap=plt.cm.Blues,title=str(dates[i]),
                showcolorbar=True, cbar_label="Water Vapor [kg/kg]")
        # plt.imshow(data[:,0,:],vmin=minval,vmax=maxval,cmap=plt.cm.Blues)
        # else:
        #     plt.imshow(data[i][:,0,:],vmin=minval,vmax=maxval+i*maxvalslope,cmap=plt.cm.Blues)
        # plt.colorbar();
        # plt.title(dates[i])
        plt.draw()
        plt.savefig("movie/qv_movie_{0:05}.png".format(i))
        # if i%100==0:
            # print(np.round(i/float(len(files))*1000)/10.0,end="% ")
            # print("{0} min={1:8.5f} max={2:8.5f}".format(dates[i],data[:,0,:].min(),data[:,0,:].max()))
            # sys.stdout.flush()
    print("done")

def main():
    """docstring for main"""
    files=glob("output/swim*")
    files.sort()
    startdate=datetime.datetime(2000,10,1,0,0)
    dates=[startdate+datetime.timedelta(i/(24.0/3)) for i in range(len(files))]
    plt_data(data=None,files=files,dates=dates)

if __name__ == '__main__':
    main()
    
