import datetime
import numpy as np

from bunch import Bunch

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
def load(filename,startyear=2004,startdate=None):
    d=np.fromfile(filename,np.float32)
    
    startlon=-122.371249999998
    nlon=2191
    dlon=0.00833333333333333
    lon=np.array([startlon+dlon*i for i in range(nlon)])

    startlat=32.9962495168051
    nlat=1291
    dlat=dlon
    lat=np.array([startlat+dlat*i for i in range(nlat)])
    
    if startdate==None:startdate=datetime.datetime(startyear,3,1,0)
    ntimes=ntimes=d.size/nlon/nlat
    dates=[startdate+datetime.timedelta(i) for i in range(ntimes)]
    
    return Bunch(data=d.reshape((ntimes,nlat,nlon)),lat=lat,lon=lon,dates=dates)
    
