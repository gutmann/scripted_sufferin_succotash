#!/usr/bin/env python
import datetime
import glob
import numpy as np
import mygis as myio
from bunch import Bunch

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
    peak_swe=max_swe(data[:180,...])
    peak_date=np.argmax(data[:180,...],axis=0) # returns first time it hits max
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


def sca(data,maxSWE=0.08):
    tmp=np.where((data>0)&(data<maxSWE))
    sca=np.zeros(data.shape)
    sca[data>0.08]=1
    if len(tmp[0])>0:
        sca[tmp]=1-(np.exp(-2.6*data[tmp]/maxSWE)-(data[tmp]/maxSWE)*np.exp(-2.6))
    return sca    

def load_noahmp(filename, startyear=2000,startdate=None):
    """Load WRF-Noahmp SWE data from a file
    """
    wrf_dir="/".join(filename.split("/")[:-1])
    geo_file=filename
    data=myio.read_nc(filename,"SNEQV").data[8*2::48,:,:]
    lat=myio.read_nc(geo_file,"lat").data
    lon=myio.read_nc(geo_file,"lon").data
    
    if startdate==None:startdate=datetime.datetime(startyear,10,1,0)
    ntimes=data.shape[0]
    
    dates=np.array([startdate+datetime.timedelta(i) for i in range(ntimes)])
    
    return Bunch(data=data,lat=lat,lon=lon,dates=dates)


def load(filename, startyear=2000,startdate=None):
    """Load WRF SWE data from a file
    """
    wrf_dir="/".join(filename.split("/")[:-1])
    geo_file=wrf_dir+"/4km_wrf_output.nc"
    data=myio.read_nc(filename,"SNOW").data
    lat=myio.read_nc(geo_file,"XLAT").data[0,...]
    lon=myio.read_nc(geo_file,"XLONG").data[0,...]
    
    if startdate==None:startdate=datetime.datetime(startyear,10,1,0)
    ntimes=data.shape[0]
    
    dates=np.array([startdate+datetime.timedelta(i) for i in range(ntimes)])
    
    return Bunch(data=data,lat=lat,lon=lon,dates=dates)


def bin_by_elevation(data,dem,mask,dz=100):
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

    return Bunch(z=z[:n]+(dz/2),veg=veg,vegmed=vegmed,vegmin=vegmin,vegmax=vegmax,
                exposed=exposed,exposedmed=exposedmed,exposedmin=exposedmin,exposedmax=exposedmax)


def load_elev_comparison(swefile="SWE_daily.nc",info="4km_wrf_output.nc",res=4,outputfile="wrf_by_elev.png",year=7,domain="FullDomain"):
    import matplotlib.pyplot as plt
    # wrf.load_elev_comparison(swefile="wrfout_d01_2008-05-01_00:00:00",res=2,outputfile="wrf_by_elev_2km_FullDomain_May2008.png")
    if res==2:
        forest=[11,12,13,14,15,18,21]
        exposed=[1,2,3,4,5,7,8,9,10,16,17,19,20,22,23,24,25,26,27]
        mayday=0
        if domain=="FullDomain":
            xmin=200
            xmax=400
            ymin=50
            ymax=400
            dz=100
        else:
            xmin=300
            xmax=380
            ymin=250
            ymax=350
            dz=200
        info=swefile
        dz=100
    else:
        forest=[1,5]
        exposed=[7,10]
        mayday=212
        for i in range(year):
            mayday+=365
        mayday+=np.floor(year/4)
        if domain=="FullDomain":
            xmin=125
            xmax=205
            ymin=80
            ymax=190
            dz=150
        else:
            xmin=175
            xmax=195
            ymin=149
            ymax=175
            dz=300
        
    print("Loading data")
    vegclass=myio.read_nc(info,"IVGTYP").data[0,...]
    mask=np.zeros(vegclass.shape)
    for f in forest:
        forested=np.where(vegclass==f)
        mask[forested]=1
    for e in exposed:
        exposed=np.where(vegclass==e)
        mask[exposed]=2
    
    dem=myio.read_nc(info,"HGT").data[0,...]
    snow=myio.read_nc(swefile,"SNOW").data[mayday,:,:]/1000.0
    
    print("Binning")
    banded=bin_by_elevation(snow[ymin:ymax,xmin:xmax],dem[ymin:ymax,xmin:xmax],mask[ymin:ymax,xmin:xmax],dz=dz)

    print("Plotting")
    plt.clf();
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
    plt.ylim(0,0.8)
    plt.ylabel("Snow Water Equivalent [m]")
    plt.xlabel("Elevation [m]")
    plt.title("WRF SWE over headwaters "+domain)
    plt.savefig(outputfile)
    
if __name__ == '__main__':
    domains=["FullDomain","GreaterFrontRange"]
    for dom in domains:
        files=glob.glob("wrfout*-05-01*")
        files.sort()
        for f in files:
            year=f.split("_")[2][:4]
            print(f,year,"2km",dom)
            load_elev_comparison(swefile=f, year=year,outputfile="wrf_by_elev_2km_MAY{}_{}.png".format(year,dom),res=2,domain=dom)
            
        for year in range(8):
            print(year,"2km",dom)
            load_elev_comparison(year=year,outputfile="wrf_by_elev_4km_MAY{}_{}.png".format(2001+year,dom),domain=dom)
    
    