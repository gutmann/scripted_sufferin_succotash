#!/usr/bin/env python
import glob
import numpy as np
from scipy.interpolate import interp1d

import mygis
from bunch import Bunch
from bilin_regrid import load_geoLUT

ncep_var="pr"
ncep_var="tas"

if ncep_var=="pr":
    ncep_files="ncep/pr/*"
    maxval=2.5
    minval=0.5
    pgw_file="ratio_file.nc"
    pgw_var="ratio"
    data_atts=Bunch(long_name="Precipitation",_FillValue=1e20, units="mm/d",GRIB_name="PRATE",
                    var_desc="Precipitation Rate",level_desc="Surface",missing_value=1e20,dataset="NMC Reanalysis")
else:
    ncep_files="ncep/tas/*.nc"
    maxval=20
    minval=-20
    pgw_file="PGW_difference_file.nc"
    pgw_var="data"
    data_atts=Bunch(long_name="2m Air Temperature",_FillValue=1e20, units="C",
                    var_desc="2m Air Temperature",level_desc="2 meters",missing_value=1e20,dataset="NMC Reanalysis")
    

month_lengths=[31,28,31,30,31,30,31,31,30,31,30,31]
months_mid_doy=np.zeros(14)
months_mid_doy[0]=-15.5
for i in range(1,len(months_mid_doy)):
    months_mid_doy[i]=months_mid_doy[i-1]+month_lengths[(i+10)%12]/2.0+month_lengths[(i+11)%12]/2.0


def load_pgw_signals():
    data=mygis.read_nc(pgw_file,pgw_var).data
    geo=mygis.read_geo(pgw_file)
    print(geo.lat.shape)
    return Bunch(data=data,geo=geo)

def load_ncep():
    """load real ncep data (by file)"""
    files=glob.glob(ncep_files)
    files.sort()
    data=[]
    times=[]
    for f in files:
        data.append(mygis.read_nc(f,ncep_var).data)
        times.append(mygis.read_nc(f,"time").data)
    
    geo=mygis.read_geo(files[0])
    return Bunch(data=data,geo=geo,files=files,time=times)

def geo_interpolate(data,geo):
    """geographic (bilinear) interpolation between the CCSM input grid and the NCEP grid"""
    geolut=load_geoLUT(data.geo.lat,data.geo.lon,geo.lat,geo.lon)
    outputdata=np.zeros((data.data.shape[0]+2,geo.lat.shape[0],geo.lat.shape[1]))
    for i in range(4):
        y=geolut[:,:,i,0].astype('i')
        x=geolut[:,:,i,1].astype('i')
        outputdata[1:-1,...]+=np.float32(data.data[:,y,x]*geolut[np.newaxis,:,:,i,2])

    outputdata[-1,...]=outputdata[1,...] #wrap around by one month (append january after december) to aid temporal interpolation
    outputdata[0,...]=outputdata[-2,...] #wrap around by one month (prepend december before january) to aid temporal interpolation
    mygis.write(ncep_var+"_pgw_test_file.nc",outputdata)
    outputdata[outputdata>maxval]=maxval
    outputdata[outputdata<minval]=minval
    return outputdata

def make_daily(pgw_signal,nceptime):
    """
    convert monthly data to daily for a given ncep year
    
    Applies a bilinear interpolation that wraps around the year 
    (e.g. Jan 1st is interpolated between Jan. and Dec. means)
    """
    f=interp1d(months_mid_doy,pgw_signal,axis=0)
    return f(nceptime)

def pgw(data,signal):
    if ncep_var=="pr":
        return data*signal
    else:
        return data+signal

def apply_pgw(ncep,pgw_signal):
    """apply pgw signal to ncep data (by file/month)"""
    ncep_pgw_signal=geo_interpolate(pgw_signal,ncep.geo) # interpolate from the ccsm pgw_signal grid to the ncep grid
    outputdata=[]
    files=[]
    ncep_daily_pgw_signal=make_daily(ncep_pgw_signal,np.arange(365)+0.5)    # interpolate from monthly means to daily values
    ncep_daily_leap_pgw_signal=make_daily(ncep_pgw_signal,365/366.*(np.arange(366)+0.5))    # same for a leap year
    for i in range(len(ncep.data)):
        if len(ncep.time[i])==365:
            outputdata.append(pgw(ncep.data[i],ncep_daily_pgw_signal))
        elif len(ncep.time[i])==366:
            outputdata.append(pgw(ncep.data[i],ncep_daily_leap_pgw_signal))
            
        files.append("pgw/"+ncep.files[i])
        
    return Bunch(data=outputdata,files=files,geo=ncep.geo,time=ncep.time)
    
def write_pgw(pgw):
    """write pgw ncep data to output files"""
    
    lat=pgw.geo.lat[:,0]
    lon=pgw.geo.lon[0,:]+360
    lat_atts =Bunch(long_name="Latitude", standard_name="latitude", units="degrees_north",
                    axis="Y",actual_range=np.array([lat.min(),lat.max()]))
    lon_atts =Bunch(long_name="Longitude", standard_name="longitude", units="degrees_east",
                    axis="X",actual_range=np.array([lon.min(),lon.max()]))
    time_atts=Bunch(long_name="Time",standard_name="time",units="hours since 1-1-1 00:00:00",
                    axis="T",avg_period="0000-00-00 06:00:00", delta_t="0000-00-00 06:00:00")
    datadims=("time","lat","lon")
    
    evars=[ Bunch(data=lat, name="lat", dtype="f",attributes=lat_atts,dims=("lat",)),
            Bunch(data=lon, name="lon", dtype="f",attributes=lon_atts,dims=("lon",)),
            Bunch(data=None,name="time",dtype="f",attributes=time_atts,dims=("time",)),
            ]
    
    for i in range(len(pgw.data)):
        times=pgw.time[i]
        evars[2].data=times
        evars[2].attributes.actual_range=np.array([times.min(),times.max()])
        
        mygis.write(pgw.files[i],pgw.data[i],dtype="f",varname=ncep_var,
                        dims=datadims,attributes=data_atts,extravars=evars)
    

def main():
    """apply monthly pgw_signals to ncep forcing data to generate a PGW NCEP"""
    print("loading pgw signals")
    pgw_signal_data=load_pgw_signals()
    
    print("loading ncep")
    ncep_data=load_ncep()
    
    print("applying PGW")
    pgw_ncep=apply_pgw(ncep_data,pgw_signal_data)
    
    print("writing data")
    write_pgw(pgw_ncep)
    

if __name__ == '__main__':
    main()