import datetime
import numpy as np
from stat_down import myio
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
