#!/usr/bin/env python
from glob import glob
import re

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from stat_down import myio as io

experiment="conus"
methodnames=dict(CAe0="BCCAr",CA="BCCA",SARe0="SAR",SDe0="BCSDd",SDmon="BCSDm",
                 CAe1="BCCAr",SARe1="SAR",SDe1="BCSDd",
                 SAR="SAR",SD="BCSDd",SDmon_c="BCSDm")

def map_vis(data,title=[],vmin=None,vmax=None,barlabel=None,cmap=None,outputdir=""):
    # print("visualizing : "+title)
    if cmap==None:
        cmap=cm.jet
    if len(data.shape)>2:
        plotdata=data[0,:,:]
    else:
        plotdata=data
    plt.clf()
    ny,nx=plotdata.shape
    geo=[35,43,-113,-101]
    if experiment=="conus":geo=[25.05,52.7,-124.7,-67.05]
    m = Basemap(projection='cyl',llcrnrlat=geo[0],urcrnrlat=geo[1],\
                llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
    
    mapimg=m.imshow(plotdata,vmin=vmin,vmax=vmax,cmap=cmap)
    if experiment=="conus":
        m.drawparallels(np.arange(25,55,5.),labels=[1,0,0,0],dashes=[1,4])
        m.drawmeridians(np.arange(-120,-65,10.),labels=[0,0,0,1],dashes=[1,4])
        m.drawstates(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawcoastlines(linewidth=0.5)
    else:
        m.drawparallels(np.arange(36,43,2.),labels=[1,0,0,0],dashes=[1,4])
        m.drawmeridians(np.arange(-112,-103,4.),labels=[0,0,0,1],dashes=[1,4])
        m.drawstates(linewidth=1.5)
    cbar=m.colorbar()
    if barlabel:
        cbar.set_label(barlabel)
    plt.title(" ".join(title))
    
    outputname=outputdir+("_".join(title)+'_map.png')
    plt.savefig(outputname)

def compare(f1,f2,dirname):
    d1=io.read_nc(f1).data
    d2=io.read_nc(f2).data
    name=f1.split("/")[1][:-3]
    delta=d1-d2
    goodpoints=delta[np.abs(delta)<1000]
    minmax=max(np.abs(goodpoints.min()),goodpoints.max())
    if re.match(".*MAP.*",name):
        deltaimg=d1-d2
        deltaimg=np.ma.array(deltaimg,mask=(deltaimg<-999)|(~np.isfinite(deltaimg)))
        map_vis(deltaimg,title=name.split("_"),cmap=cm.seismic,vmin=-450,vmax=450,outputdir="comp_"+dirname+"_pngs/")
    else:
        map_vis(d1-d2,title=name.split("_"),cmap=cm.seismic,vmin=-minmax,vmax=minmax,outputdir="comp_"+dirname+"_pngs/")

if __name__ == '__main__':
    # for test in ["sd","obs","forcing"]:
    test="sd"
    # test="obs"
    # test="forcing"
    if test=="sd":
        dir1="conusprecip_0mm/"
        dir2="conusprecip_0mm19c/"
    elif test=="obs":
        dir1="obs20c/"
        dir2="obs19c/"
    elif test=="forcing":
        dir1="forcing20c/"
        dir2="forcing19c-2/"

    d1files=glob(dir1+"*full_res_annual_MAP.nc")
    files=[]
    for f1 in d1files:
        f2=glob(dir2+f1.split("/")[1])
        if f2:
            print(f1.split("/")[1])
            try:
                compare(f1,f2[0],test)
            except Exception as e:
                print(e)
            
    