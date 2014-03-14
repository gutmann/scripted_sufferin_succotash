#!/usr/bin/env python
import glob

import numpy as np
import Nio
from bunch import Bunch

import swim_io


wrf_geo_file="/home/gutmann/data/swim/baseline/4km_wrf_output.nc"
file_search="met*.out"
output_file="noahmp_output"
variables=["SNEQV","SNOWH","FSNO"]
var_atts=Bunch(SNEQV=Bunch(units="mm",longname="Snow Water Equivalent"),
                 SNOWH=Bunch(units="m",longname="Snow Depth"),
                 FSNO=Bunch(units="[]",longname="Fractional Snow Covered Area"))
start_date="2000-10-1 00:00UTC"

def find_var(filename,varnames):
    """search through varnames until you find a variable that exists the netcdf file"""
    data=None
    for v in varnames:
        try:
            data=swim_io.read_nc(filename,v).data[0,...]
        except Exception as e:
            pass
    if data==None:
        print("Error: unable to find lat/lon variables in :"+filename)
        print("Variable names tested:")
        print(varnames)
        raise(e)
        
    return data

def load_geo(filename):
    """load lat lon data from filename"""
    latvnames=["XLAT","XLAT_M","lat","latitude"]
    lonvnames=["XLONG","XLONG_M","lon","longitude"]
    
    lat=find_var(filename,latvnames)
    lon=find_var(filename,lonvnames)
    return Bunch(lat=lat,lon=lon)

def find_pos(filename,geo):
    lat=float(filename.split("_")[-2])

    lon=float(filename.split("_")[-1][:-8])
    dists=(geo.lat-lat)**2 + (geo.lon-lon)**2
    y,x=np.unravel_index(np.argmin(dists),dists.shape)
    
    return (x,y)

def find_all_xy(files,geo):
    """find matching data points in the geo lat/lon data by reading lat,lon pos from each filename"""
    x=np.zeros(len(files))
    y=np.zeros(len(files))
    for i,f in enumerate(files):
        x[i],y[i]=find_pos(f,geo)
        
    return (x,y)

def find_times(filename):
    """return the length of the first dimension of a variable in filename"""
    d=swim_io.read_nc(filename,variables[0],returnNCvar=True)
    ntimes=d.data.shape[0]
    d.ncfile.close()
    return ntimes

def setup_nc_file(filename,lat,lon,nt,ny,nx):
    """setup a netcdf file for output"""
    
    ncfile=Nio.open_file(filename,mode="w",format="nc")
    ncfile.create_dimension("time",nt)
    ncfile.create_dimension("lat",ny)
    ncfile.create_dimension("lon",nx)

    timev=ncfile.create_variable("time","f",("time",))
    times=np.arange(nt)/48.0
    timev[:]=times.astype("f")
    timev.__setattr__("longname","Time")
    timev.__setattr__("units","Days from "+start_date)
    
    latv=ncfile.create_variable("lat","f",("lat","lon"))
    latv[:]=lat.astype("f")
    latv.__setattr__("longname","Latitude")
    latv.__setattr__("units","degrees")
    
    lonv=ncfile.create_variable("lon","f",("lat","lon"))
    lonv[:]=lon.astype("f")
    lonv.__setattr__("longname","Longitude")
    lonv.__setattr__("units","degrees")
    
    return ncfile
    
def load_data(files,varname,returndata,x,y):
    """load a grid of data, one time series form each file"""
    for i,f in enumerate(files):
        data=swim_io.read_nc(f,varname).data
        returndata[:,y[i],x[i]]=data

def write_data(data,varname,atts,ncfile):
    curvar=ncfile.create_variable(varname,"f",("time","lat","lon"))
    print(curvar.shape,curvar.typecode())
    print(data.shape,data.astype("f").dtype)
    
    curvar[:]=data.astype("f")
    try:
        print(atts)
        for k in atts.keys():
            print(k,atts[k])
            curvar.__setattr__(k,atts[k])
    except Exception as e:
        print(e)

def finalize_ncfile(ncfile):
    ncfile.close()

def main():
    geo_grid=load_geo(wrf_geo_file)
    files=glob.glob(file_search)
    files.sort()
    
    x,y=find_all_xy(files,geo_grid)

    output_lat=geo_grid.lat[y.min():y.max()+1,x.min():x.max()+1]
    output_lon=geo_grid.lon[y.min():y.max()+1,x.min():x.max()+1]
    
    x-=x.min()
    y-=y.min()
    xlen=output_lat.shape[1]
    ylen=output_lat.shape[0]
    
    ntimes=find_times(files[0])
    
    output_data=np.zeros((ntimes,ylen,xlen))
    for v in variables:
        load_data(files,v,output_data,x,y)
        
        ncfile=setup_nc_file(output_file+"_"+v,output_lat,output_lon,ntimes,ylen,xlen)
        write_data(output_data,v,var_atts[v],ncfile)
        finalize_ncfile(ncfile)
    
    

if __name__ == '__main__':
    main()