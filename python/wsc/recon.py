import datetime

import numpy as np
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
    peak_swe=max_swe(data)
    for i in range(data.shape[0]):
        nosnow=np.where(data[i,:,:]==0)
        newmelt=np.where(melt_date[nosnow]>i)
        melt_date[nosnow][newmelt]=i+1
        
        notpeak=np.where(data[i,:,:]<peak_swe)
        melt_start=np.where(peak_date[notpeak]==-1)
        peak_date[notpeak][melt_start]=i
    
    # not sure what to do if it predicts that SWE never gets to zero...
    # melt_date[melt_date==999]=data.shape[0]
    # 
    melt_time=melt_date-peak_date
    
    melt_time[melt_time<=0]=1
    
    return Bunch(rate=peak_swe/melt_time,peak=peak_swe,melt_time=melt_time,melt_date=melt_date)
        

def max_swe(data):
    """Calculate the maximum SWE at each point"""
    return data.max(axis=0)
    
    

# Assumes the following .ctl file
# DSET ^swe.dat
# UNDEF -9999
# XDEF 1950 LINEAR -112.247916666666667 0.004166666666667
# YDEF 2580 LINEAR 33.002083333333333 0.004166666666667
# ZDEF 1 LEVELS 0
# TDEF 184 LINEAR 00Z01MAR2005 1dy
# VARS 1
# SWE 0 99 SWE [m]
# ENDVARS
def load(filename,startyear=2000):
    """Load a SNODIS SWE file from disk
    
    Assumes a flat binary file as described in the comments above.
    Returns data, lat, lon, and date (based on the start year input)
    """
    d=np.fromfile(filename,np.float32)

    startlon=-112.247916666666667
    nlon=1950
    dlon=0.004166666666667
    lon=np.array([startlon+dlon*i for i in range(nlon)])

    startlat=33.002083333333
    nlat=2580
    dlat=dlon
    lat=np.array([startlat+dlat*i for i in range(nlat)])
    
    startdate=datetime.datetime(startyear,3,1,0)
    ntimes=d.size/nlon/nlat #184
    dates=[startdate+datetime.timedelta(i) for i in range(ntimes)]
    
    return Bunch(data=d.reshape((ntimes,nlat,nlon)),lat=lat,lon=lon,dates=dates)

