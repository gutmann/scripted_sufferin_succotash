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


def load(filename, startyear=2000,startdate=None,extractday=None):
    """Load WRF SWE data from a file
    """
    # find the directory in which files will be
    wrf_dir="/".join(filename.split("/")[:-1])
    if len(wrf_dir)<1:
        wrf_dir="."
        
    # geo information filename is hard coded...
    geo_file=wrf_dir+"/4km_wrf_output.nc"
    # if we only want a single day from the file, just load that day
    if extractday!=None:
        ncdata=myio.read_nc(filename,"SNOW",returnNCvar=True)
        ntimes=ncdata.data.shape[0] #store the number of times for later use
        data=ncdata.data[extractday,:,:]
        ncdata.ncfile.close()
    else:
        data=myio.read_nc(filename,"SNOW").data
        ntimes=data.shape[0]

    # load forested vs bare land cover into a mask
    vegclass=myio.read_nc(geo_file,"IVGTYP").data[0,...]
    mask=np.zeros(vegclass.shape)
    # these are the values of vegclass for forest and exposed areas in the WRF4km runs
    forest=[1,5]
    exposed=[7,10]
    # loop through LC indices setting the mask
    for f in forest:
        forested=np.where(vegclass==f)
        mask[forested]=1
    for e in exposed:
        exposed=np.where(vegclass==e)
        mask[exposed]=2
    mask[mask==0]=2 #assume unclassified areas are more exposed
    
    # load the dem and geographic data from the original file
    dem=myio.read_nc(geo_file,"HGT").data[0,...]
    lat=myio.read_nc(geo_file,"XLAT").data[0,...]
    lon=myio.read_nc(geo_file,"XLONG").data[0,...]
    
    # set up a dates array (one element if we are extracting a single day)
    if startdate==None:startdate=datetime.datetime(startyear,10,1,0)
    dates=np.array([startdate+datetime.timedelta(i) for i in range(ntimes)])
    if extractday!=None:
        dates=np.array([dates[extractday]])
    
    return Bunch(data=data/1000.0,lat=lat,lon=lon,dates=dates,dem=dem,lc=mask)

def load_base_data(swefile="SWE_daily.nc",info="4km_wrf_output.nc",res=4):
    # wrf.load_base_data(swefile="wrfout_d01_2008-05-01_00:00:00",res=2)
    
    if res==2:
        forest=[11,12,13,14,18]
        exposed=[1,2,3,4,5,7,8,9,10,16,17,19,20,22,23,24,25,26,27,15,21] #warning 15 and 21 sound like they should be forest, but on the map they are in open areas
        mayday=0
        info=swefile
        dz=100
    else:
        forest=[1,5]
        exposed=[7,10]
        mayday=212
        for i in range(year):
            mayday+=365
        mayday+=np.floor(year/4)
        
    print("Loading data")
    vegclass=myio.read_nc(info,"IVGTYP").data[0,...]
    mask=np.zeros(vegclass.shape)
    for f in forest:
        forested=np.where(vegclass==f)
        mask[forested]=1
    for e in exposed:
        exposed=np.where(vegclass==e)
        mask[exposed]=2
    lat=myio.read_nc(info,"XLAT").data[0,...]
    lon=myio.read_nc(info,"XLONG").data[0,...]
    dem=myio.read_nc(info,"HGT").data[0,...]
    snow=myio.read_nc(swefile,"SNOW").data[mayday,:,:]/1000
    return Bunch(data=snow, lat=lat,lon=lon, dem=dem, lc=mask)

def load_elev_comparison(swefile="SWE_daily.nc",info="4km_wrf_output.nc",res=4,outputfile="wrf_by_elev.png",year=7,domain="FullDomain"):
    # wrf.load_elev_comparison(swefile="wrfout_d01_2008-05-01_00:00:00",res=2,outputfile="wrf_by_elev_2km_FullDomain_May2008.png")
    import wsc.compare2lidar as c2l
    
    if res==2:
        forest=[11,12,13,14,18]
        exposed=[1,2,3,4,5,7,8,9,10,16,17,19,20,22,23,24,25,26,27,15,21] #warning 15 and 21 sound like they should be forest, but on the map they are in open areas
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
        # xmin=70;xmax=108;ymin=170;ymax=190 # Uinta Mountains
        # xmin=77;xmax=120;ymin=215;ymax=259 # Wind River Mountains
        
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
    snow=myio.read_nc(swefile,"SNOW").data[mayday,:,:]/1000
    
    print("Binning")
    banded=c2l.bin_by_elevation(snow[ymin:ymax,xmin:xmax],dem[ymin:ymax,xmin:xmax],mask[ymin:ymax,xmin:xmax],dz=dz)

    c2l.plot_elevation_bands(banded,"WRF_SWE_elev.png","WRF SWE, Elevation, and Land Cover: "+domain)
    
    
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
    
    