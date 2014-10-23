#!/usr/bin/env python
from __future__ import print_function
import glob
import os,sys
import copy
import time

import numpy as np
import nio
import netCDF4

import regrid_hi2low
from bunch import Bunch

meta=Bunch(title="WRF regridded to low-resolution projection",
    creator=os.environ['USER']+" using wrf_agg2obs.py",
    source="Model output", conventions="None")

# var_list=["Q2","RAINNC"]#,"T2","PSFC","U10","V10","SWDOWN","GLW"]
var_list=["Q2","RAINNC","T2","PSFC","U10","V10","SWDOWN","GLW"]
# samplefile="/glade/p/work/mizukami/usbr/original/downscaled/domain_UCO_12k.nc"
# samplefile="/d3/mizukami/domain_huc/domain_UCO_12k.nc"
samplefile="/glade/scratch/gutmann/usbr/domain_UCO_12k.nc"


def read_base(filename):
    """Read all variable names attibutes dimensions..."""
    f=nio.open_file(filename,mode="r")
    # f=netCDF4.Dataset(filename,mode="r")
    output=Bunch()

    variables=Bunch()
    for k in f.variables.keys():
        curvar=f.variables[k]
        data=Bunch()
        if k=="Times":
            data.data=curvar[:]
        data.atts=copy.deepcopy(curvar.__dict__)
        data.dims=copy.deepcopy(curvar.dimensions)
        data.dtype=curvar.typecode()
        variables[k]=data
        
    output.variables=variables
    output.global_attributes=copy.deepcopy(f.__dict__)
    output.dimensions=copy.deepcopy(f.dimensions)
    f.close()
    return output
    
def read_geo(filename):
    """Read geographic information (lat/lon)"""
    # f=nio.open_file(filename,mode="r")
    f=netCDF4.Dataset(filename,mode="r")
    try:
        lat=f.variables["lat"][:]
        lon=f.variables["lon"][:]
    except KeyError:
        lat=f.variables["XLAT"][0,...]
        lon=f.variables["XLONG"][0,...]
        
    lon[lon>180]-=360
    f.close()
    if len(lat.shape)==1:
        lon,lat=np.meshgrid(lon,lat)
        
    return Bunch(lat=lat,lon=lon)

def write_data(filename,data,base,geo,inputtimes=None,fillvalue=-9999):
    """Write data to filename with attributes from base and lat-lon from geo"""
    of=nio.open_file(filename,mode='w',format="nc")
    # of=netCDF4.Dataset(filename,mode='w',format="nc")
    of.title=meta.title
    of.creator=meta.creator
    of.creation_date=time.ctime()
    if base.global_attributes.has_key("history"):
        past_history=" last{"+base.global_attributes["history"]+"}"
    else:
        past_history=""
    of.history="Created : "+time.ctime()+" with wrf_agg2obs.py "+past_history
    of.source=meta.source
    of.Conventions=meta.conventions
    
    dimnames=dict(Time="Time",DateStrLen = "DateStrLen", west_east = "lon", south_north = "lat")
    longnames=dict(Q2="Specific Humidity",T2="Temperature",PSFC="Pressure",U10="Wind Speed (east-west)",
                   V10="Wind Speed (North South)",RAINNC="Precipitation",
                   SWDOWN="Downward Shortwave Radiation",GLW="Downward Longwave Radiation")
    for dimkey in base.dimensions.keys():
        if dimkey=="Time":
            of.create_dimension("Time",0)
        elif dimkey=="soil_layers_stag":
            pass
        elif dimkey=="DateStrLen":
            pass
        elif dimkey=="south_north":
            of.create_dimension(dimnames[dimkey],geo.lat.shape[0])
        elif dimkey=="west_east":
            of.create_dimension(dimnames[dimkey],geo.lat.shape[1])
        else:
            of.create_dimension(dimnames[dimkey],base.dimensions[dimkey])

    outvar=of.create_variable("Times","d",("Time",))
    outvar[:inputtimes.size]=inputtimes
    outvar=of.create_variable("lat","f",("lat","lon"))
    outvar.units="degrees"
    outvar.longname="latitude"
    outvar[:]=geo.lat.astype(np.float32)
    outvar=of.create_variable("lon","f",("lat","lon"))
    outvar[:]=geo.lon.astype(np.float32)
    outvar.units="degrees"
    outvar.longname="longitude"
    
    for var in var_list:
        curvar=base.variables[var]
        curdims=tuple([dimnames[thisdim] for thisdim in curvar.dims])
        outvar=of.create_variable(var,curvar.dtype,curdims)
        
        for k in curvar.atts.keys():
            setattr(outvar,k,curvar.atts[k])
            
        if (var=="RAINNC"):
            outvar.units="kg/m^2/s"
            outvar.description="Hourly Precipitation"
        outvar.longname=longnames[var]
        outvar._FillValue=np.float32(fillvalue)
        for i in range(data[var].shape[0]):
            outvar[i,:,:]=data[var][i,:,:].astype(curvar.dtype)
        
    of.close()
    
def read_obs(filename,var="pr"):
    """read a variable from a net cdf file"""
    f=nio.open_file(filename,mode="r")
    # f=netCDF4.Dataset(filename,mode="r")
    data=f.variables[var][:]
    f.close()
    return data

def read_wrf(filename,var="pr"):
    """read a variable from a net cdf file"""
    f=nio.open_file(filename,mode="r")
    # f=netCDF4.Dataset(filename,mode="r")
    data=f.variables[var][:]
    f.close()
    return data


def read_all_data(filelist,lastp):
    """read data from all files in filelist for all vars in var_list"""
    outputdata=Bunch()
    for f in filelist:
        print(f)
        sys.stdout.flush()
        for v in var_list:
            wrf_data=read_wrf(f,var=v)
            # convert precip from accumulated to mm/hr
            if v=="Q2":
                wrf_data=wrf_data/(1+wrf_data)
            if v=="RAINNC":
                wrf_data-=lastp
                lastp+=wrf_data[-1,:,:]
                wrf_data[1:,:,:]=np.diff(wrf_data,axis=0)
                wrf_data/=3600 #convert mm/hr to kg/s/m^2
                wrf_data[wrf_data<0]=0
            if f==filelist[0]:
                outputdata[v]=np.zeros((wrf_data.shape[0]*len(filelist),wrf_data.shape[1],wrf_data.shape[2]),dtype=np.float32)
                endpt=wrf_data.shape[0]
                outputdata[v][:endpt,:,:]=wrf_data
                startpt=endpt
            else:
                outputdata[v][startpt:endpt,:,:]=wrf_data
        startpt=endpt
        endpt=startpt+wrf_data.shape[0]
    return (outputdata,lastp)
    
def add_a_buffer(geo,buffersize):
    """Add a buffer of <buffersize> grid cells to geo.lat,geo.lon"""
    dlat=geo.lat[1,0]-geo.lat[0,0]
    dlon=geo.lon[0,1]-geo.lon[0,0]
    oldshape=geo.lat.shape
    newlat=np.zeros((oldshape[0]+buffersize*2,oldshape[1]+buffersize*2))
    newlat[buffersize:-buffersize,buffersize:-buffersize]=geo.lat
    newshape=newlat.shape
    for i in range(buffersize):
        newlat[i,:]=geo.lat[0,0]-dlat*(buffersize-i)
        newlat[newshape[0]-i-1,:]=geo.lat[-1,-1]+dlat*(buffersize-i)
    for i in range(buffersize):
        newlat[:,i]=newlat[:,buffersize]
        newlat[:,newshape[1]-i-1]=newlat[:,-buffersize-1]
    

    newlon=np.zeros((oldshape[0]+buffersize*2,oldshape[1]+buffersize*2))
    newlon[buffersize:-buffersize,buffersize:-buffersize]=geo.lon
    for i in range(buffersize):
        newlon[:,i]=geo.lon[0,0]-dlon*(buffersize-i)
        newlon[:,newshape[1]-i-1]=geo.lon[-1,-1]+dlon*(buffersize-i)
    for i in range(buffersize):
        newlon[i,:]=newlon[buffersize,:]
        newlon[newshape[0]-i-1,:]=newlon[-buffersize-1,:]
    geo.lat=newlat
    geo.lon=newlon
    
    

def main(root_dir=None,var=None,filesearch=None):
    # WRF hourly 4km files should be in : 
    # /glade/p/ral/RHAP/asd000/HW2010.2/
    #      pgw/wrfout/wrfout*
    #      ctrl/wrfout/wrfout*
    wrffilesearch="*.nc"
    if var==None:
        var=var_list
    if filesearch:
        wrffilesearch=filesearch
    if root_dir:
        wrffilesearch=root_dir+wrffilesearch
    else:
        root_dir=""
#   wrffilesearch="*.nc"
    buffersize=2
    geo=read_geo(samplefile)
    files=glob.glob(wrffilesearch)
    files.sort()
    wrfgeo=read_geo(files[0])
    # print("reading WRF base data")
    wrf_base=read_base(files[0])
    
    if buffersize>0:
        add_a_buffer(geo,buffersize)
    
    lon,lat=wrfgeo.lon,wrfgeo.lat
    print("computing Geo LUT")
    geoLUT=regrid_hi2low.load_geoLUT(lat,lon,geo.lat,geo.lon)
    if buffersize>0:
        geo.lat=geo.lat[buffersize:-buffersize,buffersize:-buffersize]
        geo.lon=geo.lon[buffersize:-buffersize,buffersize:-buffersize]
    fillvalue=-9999 #base.variables.pr.atts["_FillValue"][0]
    lastp=np.zeros(wrfgeo.lon.shape)
    years=range(2000,2009)
    starttime=0.0
    months=["01","02","03","04","05","06","07","08","09","10","11","12"]
    for y in years:
        for m in months:
            curfiles=glob.glob("*"+str(y)+"-"+m+"*.nc")
            if len(curfiles)>0:
                curfiles.sort()
                outputdata,lastp=read_all_data(curfiles,lastp)
                outputtimes=np.arange(outputdata[var_list[0]].shape[0])/24.0+starttime
                starttime=outputtimes[-1]+1.0/24.0
                for v in var_list:
                    print(v,end=" ")
                    outputdata[v]=regrid_hi2low.regrid_hi2low(outputdata[v],geoLUT=geoLUT,FillValue=fillvalue).data
                    if buffersize>0:
                        outputdata[v]=outputdata[v][:,buffersize:-buffersize,buffersize:-buffersize]
                outputdir="/glade/u/home/gutmann/scratch/usbr/wrf_aggregation/"
                outputfile=outputdir+"regridded_wrf_output_"+str(y)+"_"+m
                print("  writing...")
                write_data(root_dir+outputfile,outputdata,wrf_base,geo,outputtimes)

if __name__ == '__main__':
    main()
