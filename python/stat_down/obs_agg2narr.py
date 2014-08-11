#!/usr/bin/env python
import glob
import os
import copy
import time
import numpy as np
import mygis
import regrid_hi2low
from bunch import Bunch

# process_var="tasmin"
process_var="tasmax"
# process_var="tas"
# process_var="pr"

# narr_samplefile="/Volumes/G-SAFE/usbr/statistical/narr_sample.nc"
def get_sample_file(process_var):
    narr_samplefile=None
    if process_var=="pr":
        narr_samplefile="/d5/gutmann/cc-downscaling-test/ccsm/pr/prate.run5.complete.1960-2080.nc"
    elif process_var=="tas":
        narr_samplefile="/d5/gutmann/cc-downscaling-test/ccsm/tas/tas.run5.complete.1960-2080.nc"
    elif process_var=="tasmin":
        narr_samplefile="/d5/gutmann/cc-downscaling-test/ccsm/tasmin/tmin.run5.complete.1900-2080.nc"
    elif process_var=="tasmax":
        narr_samplefile="/d5/gutmann/cc-downscaling-test/ccsm/tasmax/tmax.run5.complete.1900-2080.nc"
    return narr_samplefile
# obs_filesearch="*"+process_var+".*.nc"
meta=Bunch(title="Observed regridded to model forcing projection",
    creator=os.environ['USER']+" using obs_agg2narr.py",
    source="UW via USBR", conventions="None")

# obs_filesearch="nldas_met*.nc"
# meta=Bunch(title="Observed regridded to model forcing projection",
    # creator=os.environ['USER']+" using obs_agg2narr.py",
    # source="Maurer via USBR", conventions="None")


def read_base(filename):
    f=mygis.Dataset(filename,mode="r")
    output=Bunch()

    variables=Bunch()
    for k in f.variables.keys():
        curvar=f.variables[k]
        if (k!="pr") and (k!="tasmin") and (k!="tasmax") and (k!="TREFMN") and (k!="TREFMX") and (k!="tas"):
            data=Bunch(data=curvar[:])
        else:
            data=Bunch()
        data.atts=copy.deepcopy(curvar.__dict__)
        # data.dims=copy.deepcopy(curvar.dimensions)
        data.dims=curvar.dimensions
        data.dtype=curvar.dtype
        variables[k]=data
        
    output.variables=variables
    output.global_attributes=copy.deepcopy(f.__dict__)
    # output.dimensions=copy.deepcopy(f.dimensions)
    output.dimensions=f.dimensions
    output.ncfile=f
    # f.close()
    return output

def write_data(filename,data,nbase,obase,outputvar="pr",fillvalue=None):
    of=mygis.Dataset(filename,mode='w')
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
            of.createDimension("time",0)
        else:
            of.createDimension(dimkey,int(str(nbase.dimensions[dimkey]).split("=")[-1]))
            # of.createDimension(dimkey,int(nbase.dimensions[dimkey]))
    for var in nbase.variables.keys():
        
        # Create the variable in the output netcdf file
        curvar=nbase.variables[var]
        ovar=var
        if var=="TREFMN":ovar="tasmin"
        if var=="TREFMX":ovar="tasmax"
        if (var!="pr") and (var!="tasmin") and (var!="tasmax") and (var!="TREFMN") and (var!="TREFMX") and (var!="tas"):
            if var=="time":
                try:
                    curvar=obase.variables[var]
                    outvar=of.createVariable(var,curvar.dtype,tuple(curvar.dims))
                except:
                    print("ERROR managing time")
            elif var=="time_bnds":
                pass
            else:
                outvar=of.createVariable(ovar,curvar.dtype,tuple(curvar.dims))
        else:
            curvar=obase.variables[ovar]
            outputdims=[]
            for thisdim in curvar.dims:
                if thisdim=="latitude":
                    thisdim="lat"
                if thisdim=="longitude":
                    thisdim="lon"
                outputdims.append(thisdim)
                    
            outvar=of.createVariable(outputvar,curvar.dtype,tuple(outputdims),fill_value=fillvalue)
            
        # copy attributes in
        if (var!="time_bnds"):
            for k in curvar.atts.keys():
                if k!="_FillValue":
                    setattr(outvar,k,curvar.atts[k])
        
        # finally copy the data in
        if (var=="pr") or (var=="tasmin") or (var=="tasmax") or (var=="TREFMN") or (var=="TREFMX") or (var=="tas"):
            try:
                setattr(outvar,"units",getattr(obase.variables,outputvar).atts["units"])
            except:
                pass
            for i in range(data.shape[0]):
                outvar[i,:,:]=data[i,:,:].astype(curvar.dtype)
        elif var=="time":
            obs_time=obase.variables.time.data
            for i in range(obs_time.size):
                outvar[i]=obs_time[i].astype(curvar.dtype)
        elif var=="time_bnds":
            pass
        else:
            outvar[:]=curvar.data[:].astype(curvar.dtype)
    of.close()
    
def read_obs(filename,var="pr"):
    return mygis.read_nc(filename,var).data
                
def main(root_dir=None,var=None,filesearch=None):
    if var==None:
        var=process_var
    obs_filesearch="*"+var+".*.nc"
    if filesearch:
        obs_filesearch=filesearch
    if root_dir:
        obs_filesearch=root_dir+obs_filesearch
    else:
        root_dir=""
    
    files=glob.glob(obs_filesearch)
    obsgeo=mygis.read_geo(files[0])
    if len(obsgeo.lon.shape)==1:
        obslon,obslat=np.meshgrid(obsgeo.lon,obsgeo.lat)
    else:
        obslon,obslat=(obsgeo.lon,obsgeo.lat)
    
    narr_samplefile=get_sample_file(var)
    narrgeo=mygis.read_geo(narr_samplefile)
    narrbase=read_base(narr_samplefile)
    print("computing Geo LUT")
    geoLUT=regrid_hi2low.load_geoLUT(obslat,obslon,narrgeo.lat,narrgeo.lon)
        
    try:
        fillvalue=narrbase.variables.pr.atts["_FillValue"][0]
    except:
        fillvalue=1e20
    for f in files:
        print("Working on : "+f)
        obs_data=read_obs(f,var=var)
        obs_base=read_base(f)
        outputdata=regrid_hi2low.regrid_hi2low(obs_data,geoLUT=geoLUT,FillValue=fillvalue)
        outputfile="regridded_"+f.split('/')[-1]
        print("  writing...")
        write_data(root_dir+outputfile,outputdata.data,narrbase,obs_base,outputvar=var,fillvalue=fillvalue)

if __name__ == '__main__':
    main()