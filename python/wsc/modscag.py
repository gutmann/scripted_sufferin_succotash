import datetime
import date_fun
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
    # find the maximum SCA in the first half of the year (so a peak in Dec doesn't screw up melt "rates")
    peak_sca=max_sca(data[:data.shape[0]/2,:,:])
    # loop through time finding times that have reached "peak" or that have melted
    # record first melt date and last peak date
    for i in range(data.shape[0]):
        nosnow=np.where(data[i,:,:]==0)
        newmelt=np.where(melt_date[nosnow]>i)
        melt_date[nosnow][newmelt]=i+1
        
        notpeak=np.where(data[i,:,:]<peak_sca)
        melt_start=np.where(peak_date[notpeak]==-1)
        peak_date[notpeak][melt_start]=i
    
    # not sure what to do if it predicts that SCA never gets to zero...
    melt_date[melt_date==999]=data.shape[0]
    
    melt_time=melt_date-peak_date
    # if meltdate is before peak date just set it to 1
    melt_time[melt_time<=0]=1
    
    return Bunch(peak=peak_sca,melt_time=melt_time,melt_date=melt_date)
        

def max_sca(data):
    """Calculate the maximum SCA at each point (quite likely 1?)"""
    return data.max(axis=0)
    
def fill(data):
    """Fill missing MODSCAG values by linear interpolation in time. 
    
    Some data in MODSCAG are missing (-9999 or >100 (various flags, water, cloud, ...))
    """
    def bad_data(datavalue):
        return (datavalue<0) or (datavalue>1)
    
    # loop through time
    for i in range(data.shape[0]):
        # find all missing values in the current time step
        tmp=np.where((data[i,:,:]<0)|(data[i,:,:]>1))
        # loop through missing values
        for j in range(len(tmp[0])):
            # current spatial position to test
            x=tmp[0][j]
            y=tmp[1][j]
            # backwards search index starting point
            test_i=i-1
            # note, after i=0 test_i should always = i-1
            if (i==0):
                # search backwards for a good datapoint
                while bad_data(data[test_i,x,y]):test_i-=1
            i_start=test_i
            # forward search index starting point
            test_i=i+1
            # search forwards for a good datapoint
            while bad_data(data[test_i,x,y]):test_i+=1
            
            # linear interpolation between last good point and next good point
            linterp=(float(test_i)-i)/(float(test_i)-i_start)
            data[i,x,y]=data[test_i,x,y]*(1-linterp)+data[i_start,x,y]*linterp
    
    # return data
    
def month_str2num(month):
    """convert a 3char month name to the number
    
    e.g. JAN = 1, FEB = 2, etc.
    """
    datenumbers=dict(JAN=1,FEB=2,MAR=3,APR= 4,MAY= 5,JUN= 6,
                     JUL=7,AUG=8,SEP=9,OCT=10,NOV=11,DEC=12)
    return datenumbers[month.upper()]

def read_good_days(filename="dateselect.txt"):
    """Read dates QCed as good from filename (MODSCAGdir/dateselect.txt)
    
    dateselect.txt has one line per date in the format """
    with open(filename) as f:
        dates=[]
        indices=[]
        for l in f:
            day=int(l[0:2])
            month=month_str2num(l[2:5])
            year=int(l[5:])
            dates.append(datetime.datetime(year,month,day))
            indices.append(int(np.round(date_fun.date2mjd(year,month,day,0,0)
                                       -date_fun.date2mjd(year,1,1,0,0))))
    return Bunch(dates=dates,indices=indices)
            
# Assumes the following .ctl file
# DSET ^fsca.dat
# UNDEF -9999
# XDEF 1950 LINEAR -112.247916666666667 0.004166666666667
# YDEF 2580 LINEAR 33.002083333333333 0.004166666666667
# ZDEF 1 LEVELS 0
# TDEF 366 LINEAR 00Z01JAN200X 1dy
# VARS 1
# SCA 0 1 SCA []
# ENDVARS
def load(filename,startyear=2008,all_days=False):
    """Load a MODSCAG SCA file from disk
    
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
    
    startdate=datetime.datetime(startyear,1,1,0)
    ntimes=d.size/nlon/nlat #366
    dates=[startdate+datetime.timedelta(i) for i in range(ntimes)]
    
    dateselectfile="/".join(filename.split("/")[:-1])+"/dateselect.txt"
    if all_days:
        gooddates=np.arange(ntimes)
    else:
        gooddates=read_good_days(dateselectfile)
    
    # fill(data)
    return Bunch(data=d.reshape((ntimes,nlat,nlon)),lat=lat,lon=lon,dates=dates,gooddates=gooddates)

