#!/usr/bin/env python
import glob,os,time,sys

import numpy as np
import netCDF4

from bunch import Bunch

def get_data(filename):
    """Load temperature, time, lat,lon data and attributes into a meaningful structure"""
    dataset=netCDF4.Dataset(filename,"r")
    data=dataset.variables["TREFHT"][:]
    return Bunch(ncdata=dataset,tref=data)

def compute_daily_func(data,func=np.min):
    sz=data.tref.shape
    outputdata=func(data.tref.reshape((sz[0]/4,4,sz[1],sz[2])), axis=1)
    return Bunch(ncdata=data.ncdata,data=outputdata)

def copy_atts(inputdata,output,addhistory=""):
    for att in inputdata.ncattrs():
        if att=="history":
            output.setncattr(att,addhistory+inputdata.getncattr(att))
        else:
            output.setncattr(att,inputdata.getncattr(att))

def copy_var(ncinput_var,ncoutput_var,step=1):
    """docstring for copy_var"""
    ncoutput_var[:]=ncinput_var[::step]
    copy_atts(ncinput_var,ncoutput_var)
    

def write_data(filename,data):
    nc=netCDF4.Dataset(filename,mode="w",clobber=True,format="NETCDF4")
    nc.createDimension("time",None)
    nc.createDimension("lat",data.data.shape[1])
    nc.createDimension("lon",data.data.shape[2])
    
    times=nc.createVariable("time","d",("time",))
    lats=nc.createVariable("lat","d",("lat",))
    lons=nc.createVariable("lon","d",("lon",))
    odata=nc.createVariable(data.name,"f",("time","lat","lon"))
    
    copy_var(data.ncdata.variables["time"],times,step=4)
    copy_var(data.ncdata.variables["lat"],lats)
    copy_var(data.ncdata.variables["lon"],lons)
    
    odata[:]=data.data
    odata.setncatts(data.attributes)
    
    newhistory="Daily data created:"+time.strftime("%Y/%m/%d %H:%M:%S")+\
                " by:"+os.environ["USER"]+\
                " using script:"+__file__+"\n"
    copy_atts(data.ncdata,nc,addhistory=newhistory)
    nc.close()
    
def print_name():
    print(__name__)
    

def main():
    files=glob.glob("b30*TREFHT*")
    files.sort()
    for f in files:
        print(f)
        data=get_data(f)
        
        tmin=compute_daily_func(data,np.min)
        tmin.name="TREFMN"
        tmin.attributes=dict(long_name="Minimum daily reference height temperature",
                             units="K")
        outputfile=f.replace("TREFHT","TREFMN")
        write_data(outputfile,tmin)

        tmax=compute_daily_func(data,np.max)
        tmax.name="TREFMX"
        tmax.attributes=dict(long_name="Maximum daily reference height temperature",
                             units="K")
        outputfile=f.replace("TREFHT","TREFMX")
        write_data(outputfile,tmax)
        
        data.ncdata.close()


if __name__ == '__main__':
    main()