import numpy as np
import datetime
from bunch import Bunch

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

