#!/usr/bin/env python
import glob

import numpy as np

# from stat_down import myio
import mygis as myio
import date_fun
from bunch import Bunch

# CHANGE HERE (and below in main)
search_dir="SAR3conus/"
# END CHANGE HERE


stat_data_dir="/d2/gutmann/usbr/stat_data/DAILY/down/"
maskfile_6km=stat_data_dir+"CA/ncep/tasmax/BCCA_6km_T62_tmax.gauss.2m.2000.nc"
maskfile_12km=stat_data_dir+"CA/ncep/tasmax/BCCA_12km_T62_tmax.gauss.2m.2000.nc"

dims=("time","lat","lon")
FILL_VALUE=np.float32(-9999) # if this doesn't work, (np.zeros(3)-9999).astype("f")[1] should work...
lat_atts=dict(units="degrees_north",axis="Y", long_name="Latitude",standard_name="latitude")
lat_info=Bunch(data=None,name="lat",dims=("lat",),dtype="f",attributes=lat_atts)

lon_atts=dict(units="degrees_east",axis="X", long_name="Longitude",standard_name="longitude")
lon_info=Bunch(data=None,name="lon",dims=("lon",),dtype="f",attributes=lon_atts)

time_atts=dict(units="days since 1940-01-01 00:00:00",axis="T",long_name="Time",standard_name="time",calendar="standard")
time_info=Bunch(data=None,name="time",dims=("time",),dtype="d",attributes=time_atts)

ta_atts=dict(units="C",_FillValue=FILL_VALUE,missing_value=FILL_VALUE)
tas_info=Bunch(data=None,name="tas",dims=dims,dtype="f",attributes=ta_atts)
tasmax_info=Bunch(data=None,name="tasmax",dims=dims,dtype="f",attributes=ta_atts)
tasmin_info=Bunch(data=None,name="tasmin",dims=dims,dtype="f",attributes=ta_atts)

pr_atts=dict(units="mm/day",_FillValue=FILL_VALUE,missing_value=FILL_VALUE)
pr_info=Bunch(data=None,name="pr",dims=dims,dtype="f",attributes=pr_atts)

data_info=dict(pr=pr_info,tasmax=tasmax_info,tasmin=tasmin_info,tas=tas_info)


def load_mask(res="12km"):
    maskfile=maskfile_12km if res=="12km" else maskfile_6km
    return myio.read_nc(maskfile,"tasmax").data.mask[0,...]

def load_geo(res="12km"):
    geofile=maskfile_12km if res=="12km" else maskfile_6km
    lat=myio.read_nc(geofile,"lat").data
    lon=myio.read_nc(geofile,"lon").data
    return lat,lon

def time_gen(year,model):
    t_base=date_fun.date2mjd(1940,1,1,0,0,0)
    t_start=date_fun.date2mjd(year,1,1,0,0,0)-t_base
    if model=="ccsm":
        t_stop=t_start+365
    else:
        t_stop=date_fun.date2mjd(year+1,1,1,0,0,0)-t_base
    return np.arange(t_start,t_stop)

def update_files_for_year(year,res,variable,model,mask):
    files=glob.glob(search_dir+model+"/"+variable+"/BCSAR*"+res+"*"+str(year)+"*")
    files.sort()
    time_info.data=time_gen(year,model)
    data=np.concatenate(myio.read_files(files,variable))
    for i in range(data.shape[0]):
        data[i,...][mask]=FILL_VALUE
    
    info=data_info[variable]
    newfilename=files[0].replace(".nc","")[:-3]+".nc"
    print(newfilename)
    extra_vars=[lat_info,lon_info,time_info]
    myio.write(newfilename,data,varname=info.name,dtype=info.dtype,dims=info.dims,
               attributes=info.attributes,extravars=extra_vars)

def main():
    # res=["12km","6km"]
    # variables=["pr","tasmax","tasmin"]
    # models=["ncep","narr","ccsm"]
    
    # CHANGE HERE
    res=["6km"]
    models=["ccsm"]
    # variables=["tasmax","tasmin"]
    variables=["tas","pr"]
    # END CHANGE HERE
    
    
    for r in res:
        mask=load_mask(r)
        lat,lon=load_geo(r)
        lat_info.data=lat
        lon_info.data=lon
        for m in models:
            if m=="ccsm":
                years=range(1979,2060)
            else:
                years=range(1979,2009)
            for v in variables:
                for y in years:
                    print(r,m,v,y)
                    # try:
                    update_files_for_year(year=y,res=r,variable=v,model=m,mask=mask)
                    # except Exception as e:
                    #     print(e)
                    #     print(r,m,v,y)

if __name__ == '__main__':
    main()