#!/usr/bin/env python
import datetime

import numpy as np
from bunch import Bunch

import mygis as myio

FILLVALUE=-9999

def stats(data):
    """Calculate the rate of melt from peak to 0
    
    Assumes that data starts at peak and decreases from there
    Takes the first point data dips below peak as the onset of melt
    Takes the first day data get to 0 as meltout
    
    Rate is calculated as peak/time [m/day]
    Also returns peak swe [m], melt time [days], melt DOY [days]
    """
    melt_date=np.zeros(data.shape[1:])+999
    peak_date=np.zeros(data.shape[1:])-1
    peak_swe=max_swe(data)
    peak_date=np.argmax(data,axis=0) # returns first time it hits max
    melt_date=np.argmin(data,axis=0) # returns first time it hits min (0)
    
    # not sure what to do if it predicts that SWE never gets to zero...
    # melt_date[melt_date==data.shape[0]]=data.shape[0]
    # 
    melt_time=melt_date-peak_date
    
    melt_time[melt_time<=0]=1
    
    return Bunch(rate=peak_swe/melt_time,peak=peak_swe,melt_time=melt_time,melt_date=melt_date)
    

def max_swe(data):
    """Calculate the maximum SWE at each point"""
    return data.max(axis=0)

def fill_first(data):
    """Fill in the first time slice in data with the first non-missing value"""
    tmp=np.where(data[0,...]==FILLVALUE)
    for y,x in zip(tmp[0],tmp[1]):
        good_values=np.where(data[:,y,x]!=FILLVALUE)
        if len(good_values[0])>0:
            data[0,y,x]=data[good_values[0][0],y,x]
    

def fill_missing(data):
    """Fill in missing (-9999) data using previous value"""
    
    print("Filling in missing values")
    fill_first(data)
    for i in range(1,data.shape[0]):
        tmp=np.where(data[i,...]==FILLVALUE)
        if len(tmp[0])>0:
            data[i,...][tmp]=data[i-1,...][tmp]



# Assumes the following .ctl file
# DSET ^%y4/SWE/SWE_Daily0600UTC_WesternUS_%y4.dat
# UNDEF -9999
# options template
# XDEF 2191 LINEAR -122.371249999998 0.00833333333333333
# YDEF 1291 LINEAR 32.9962495168051 0.00833333333333333
# ZDEF 1 LEVELS 1000
# tdef 3259 linear 00z30sep2003 1dy
# VARS 1
# swe 0 99 Snow Water Equivalant (Unit: meter); Snapshot at 0600UTC; Arbitrarily assigned to 0000UTC of the same day.
# ENDVARS
def load(filename,startyear=2004,startdate=None,fill=True):
    """Load a SNODAS SWE file
    
    Assumes a flat binary file as described above, but calculates the number of days present
    Returns the data along with lat, lon, and date (based on startdate or startyear)
    """
    if filename!=None:
        d=np.fromfile(filename,np.float32)
    else:
        d=None
    
    startlon=-122.371249999998
    nlon=2191
    dlon=0.00833333333333333
    lon=np.array([startlon+dlon*i for i in range(nlon)])

    startlat=32.9962495168051
    nlat=1291
    dlat=dlon
    lat=np.array([startlat+dlat*i for i in range(nlat)])
    
    if startdate==None:startdate=datetime.datetime(startyear,1,1,0)
    if d==None:
        ntimes=365
    else:
        ntimes=ntimes=d.size/nlon/nlat
    dates=[startdate+datetime.timedelta(i) for i in range(ntimes)]
    
    if d!=None:
        d=d.reshape((ntimes,nlat,nlon))
        if fill:
            fill_missing(d)
    
    return Bunch(data=d,lat=lat,lon=lon,dates=dates)
    

def bin_by_elevation(data,dem,mask):
    """docstring for bin_by_elevation"""
    minz=dem[dem>100].min()
    maxz=dem[dem<5000].max()
    dz=100
    
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
        curexp=np.where((dem>z[i])&(dem<=z[i+1])&(mask==2)&np.isfinite(data)&(data>0)&(data<20))
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


def load_elev_comparison(swefile="SWE_Daily0600UTC_WesternUS_2010.dat"):
    import matplotlib.pyplot as plt
    
    demfile="snodas_dem.nc"
    forestfile="snodas_forest.nc"
    forest=[1]
    bare=[0]
    
    mayday=120
    minx=1800; miny=600; maxx=None;maxy=1100
    
    print("Loading data")
    vegclass=myio.read_nc(forestfile).data[miny:maxy,minx:maxx]
    forested=np.where(vegclass==forest[0])
    exposed=np.where(vegclass==bare[0])
    mask=np.zeros(vegclass.shape)
    mask[forested]=1
    mask[exposed]=2
    
    dem=myio.read_nc(demfile).data[miny:maxy,minx:maxx]
    snodas=load(swefile)
    snow=snodas.data[mayday,miny:maxy,minx:maxx]
    
    print("Binning")
    banded=bin_by_elevation(snow,dem,mask)

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
    plt.ylim(0,0.7)
    plt.ylabel("Snow Water Equivalent [m]")
    plt.xlabel("Elevation [m]")
    plt.title("SNODAS SWE over front range+")
    plt.savefig("snodas_by_elev.png")
    
if __name__ == '__main__':
    load_elev_comparison()
    
    