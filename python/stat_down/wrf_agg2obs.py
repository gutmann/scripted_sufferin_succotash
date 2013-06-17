#!/usr/bin/env python
import glob
import os
import copy
import time
import numpy as np
import nio
import regrid_hi2low
from bunch import Bunch

# process_var="tasmin"
process_var="tasmax"
# process_var="pr"

meta=Bunch(title="WRF regridded to low-resolution projection",
    creator=os.environ['USER']+" using wrf_agg2narr.py",
    source="Model output", conventions="None")



def read_base(filename):
    f=nio.open_file(filename,mode="r")
    output=Bunch()

    variables=Bunch()
    for k in f.variables.keys():
        curvar=f.variables[k]
        if (k!="pr") and (k!="tasmin") and (k!="tasmax"):
            data=Bunch(data=curvar[:])
        else:
            data=Bunch()
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
    f=nio.open_file(filename,mode="r")
    lat=f.variables["lat"][:]
    lon=f.variables["lon"][:]
    if lon.max()>180:
        lon=lon-360
    f.close()
    return Bunch(lat=lat,lon=lon)

def write_data(filename,data,nbase,obase,outputvar="pr"):
    of=nio.open_file(filename,mode='w')
    of.title=meta.title
    of.creator=meta.creator
    of.creation_date=time.ctime()
    if obase.global_attributes.has_key("history"):
        past_history=" last{"+obase.global_attributes["history"]+"}"
    else:
        past_history=""
    of.history="Created : "+time.ctime()+past_history
    of.source=meta.source
    of.Conventions=meta.conventions
    
    for dimkey in nbase.dimensions.keys():
        if dimkey=="time":
            of.create_dimension("time",0)
        else:
            of.create_dimension(dimkey,nbase.dimensions[dimkey])
    
    for var in nbase.variables.keys():
        curvar=nbase.variables[var]
        if (var!="pr"):
            outvar=of.create_variable(var,curvar.dtype,tuple(curvar.dims))
        else:
            outvar=of.create_variable(outputvar,curvar.dtype,tuple(curvar.dims))
        for k in curvar.atts.keys():
            setattr(outvar,k,curvar.atts[k])
        if (var=="pr"): # or (var=="tasmin") or (var=="tasmax"):
            setattr(outvar,"units",getattr(obase.variables,outputvar).atts["units"])
            for i in range(data.shape[0]):
                outvar[i,:,:]=data[i,:,:].astype(curvar.dtype)
        elif var=="time":
            obs_time=obase.variables.time.data
            for i in range(obs_time.size):
                outvar[i]=obs_time[i].astype(curvar.dtype)
        else:
            outvar[:]=curvar.data[:].astype(curvar.dtype)
    of.close()
    
def read_obs(filename,var="pr"):
    f=nio.open_file(filename,mode="r")
    data=f.variables[var][:]
    f.close()
    return data
            
def main(root_dir=None,var=None,filesearch=None):
    wrffilesearch="*"+process_var+".*.nc"
    samplefile="nldas_met*.nc"
    if var==None:
        var=process_var
    if filesearch:
        wrffilesearch=filesearch
    if root_dir:
        wrffilesearch=root_dir+wrffilesearch
    else:
        root_dir=""
    lowgeo=read_geo(samplefile)
    lowbase=read_base(samplefile)
    files=glob.glob(wrffilesearch)
    geo=read_geo(files[0])
    lon,lat=np.meshgrid(geo.lon,geo.lat)
    print("computing Geo LUT")
    geoLUT=regrid_hi2low.load_geoLUT(lat,lon,geo.lat,geo.lon)
    fillvalue=base.variables.pr.atts["_FillValue"][0]
    for f in files:
        print("Working on : "+f)
        wrf_data=read_wrf(f,var=var)
        wrf_base=read_base(f)
        outputdata=regrid_hi2low.regrid_hi2low(wrf_data,geoLUT=geoLUT,FillValue=fillvalue)
        outputfile="regridded_"+f.split('/')[-1]
        print("  writing...")
        write_data(root_dir+outputfile,outputdata.data,base,wrf_base,outputvar=var)

if __name__ == '__main__':
    main()
