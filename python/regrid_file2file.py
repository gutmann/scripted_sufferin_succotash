#!/usr/bin/env python
import os,glob

import numpy as np

import mygis
from bunch import Bunch
from bilin_regrid import regrid,load_geoLUT

def main(f1=None,f2=None,output_dir="./",inputvar="tos",outputvar="sst"):
    """docstring for main"""
    if f1!=None:
        inputfile=f1
    else:
        inputfile="tos_day_ACCESS1-3_historical_r1i1p1_19500101-19591231.nc"
    if f2!=None:
        reference=f2
    else:
        reference="../../completed/historical/access/subset/ta_6hrLev_ACCESS1-3_historical_r1i1p1_1950010106-1951010100.nc"
    
    files=glob.glob(inputfile)
    files.sort()
    
    print("Reading geo data")
    geo1=mygis.read_geo(files[0])
    geo2=mygis.read_geo(reference)
    
    latatts=Bunch(standard_name="latitude",long_name="latitude coordinate",units="degrees_north")
    lonatts=Bunch(standard_name="longitude",long_name="longitude coordinate",units="degrees_east")
    timeatts=None
    extra_vars=[Bunch(data=geo2.lat,name="lat",dims=("lat","lon"),dtype="f",attributes=latatts),
                Bunch(data=geo2.lon,name="lon",dims=("lat","lon"),dtype="f",attributes=lonatts),
                Bunch(data=None,name="time",dims=("time",),dtype="d",attributes=timeatts)
                ]
    
    missing_value=0
    ncdata=mygis.read_nc(files[0],inputvar,returnNCvar=True)
    geoLUT=load_geoLUT(lat1=geo1.lat,lon1=geo1.lon,lat2=geo2.lat,lon2=geo2.lon,
                        mask=(ncdata.data[0,...]==missing_value),winhalfsize=7)
    ncdata.ncfile.close()
    for f in files:
        print("Reading data: "+f)
        data=mygis.read_nc(f,inputvar).data
        data_atts=mygis.read_atts(f,inputvar)
        global_atts=mygis.read_atts(f,global_atts=True)
        time_atts=mygis.read_atts(f,"time")
        
        print("Regridding")
        output=regrid(data,geoLUT=geoLUT,missing=0)
    
        print("Writing data")
        extra_vars[-1].data=mygis.read_nc(f,"time").data
        extra_vars[-1].attributes=time_atts
        
        mygis.write(output_dir+f.replace(inputvar,outputvar),output,dims=("time","lat","lon"),varname=outputvar,attributes=data_atts,
                    global_attributes=global_atts,extravars=extra_vars)
    
if __name__ == '__main__':
    
    
    if not os.path.exists("subset"):
        os.mkdir("subset")
        
    subsetlist=[
            "../../completed/historical/access/subset/ta_6hrLev_ACCESS1-3_historical_r1i1p1_1950010106-1951010100.nc",
            "../../completed/historical/bcc/subset/ta_6hrLev_bcc-csm1-1-m_historical_r1i1p1_195001010000-195003141800.nc",
            "../../completed/historical/bnu/subset/ta_6hrLev_BNU-ESM_historical_r1i1p1_1950010100-2005123118.nc",
            "../../completed/historical/canesm/subset/ta_6hrLev_CanESM2_historical_r1i1p1_195001010000-195012311800.nc",
            "../../completed/historical/ccsm/subset/ta_6hrLev_CCSM4_historical_r6i1p1_1950010106-1950033118.nc",
            "../../completed/historical/cnrm_cm5/subset/ta_6hrLev_CNRM-CM5_historical_r1i1p1_195001010600-195002010000.nc",
            "../../completed/historical/fgoals/subset/ta_6hrLev_FGOALS-g2_historical_r1i1p1_1950010106-1951010100.nc",
            "../../completed/historical/gfdl_cm3/subset/ta_6hrLev_GFDL-CM3_historical_r1i1p1_1950010100-1950123123.nc",
            "../../completed/historical/gfdl_esm/subset/ta_6hrLev_GFDL-ESM2M_historical_r1i1p1_1951010100-1955123123.nc",
            "../../completed/historical/ipsl_mr/subset/ta_6hrLev_IPSL-CM5A-MR_historical_r1i1p1_1950010103-1959123121.nc",
            "../../completed/historical/giss_e2h/subset/ta_6hrLev_GISS-E2-H_historical_r6i1p3_195001010600-195007010000.nc",
            "../../completed/historical/miroc5/subset/ta_6hrLev_MIROC5_historical_r1i1p1_1950010100-1950013118.nc",
            "../../completed/historical/miroc_esm/subset/ta_6hrLev_MIROC-ESM_historical_r1i1p1_1950010106-1950020100.nc",
            "../../completed/historical/mk3/subset/ta_6hrLev_CSIRO-Mk3-6-0_historical_r1i1p1_195001010600-195101010000.nc",
            "../../completed/historical/mri_cgcm3/subset/ta_6hrLev_MRI-CGCM3_historical_r1i1p1_195001010000-195001311800.nc",
            "../../completed/historical/noresm/subset/ta_6hrLev_NorESM1-M_historical_r1i1p1_1950010100-1950063018.nc"]
        
    toslist=[
        "tos_day_ACCESS1-3_historical_r1i1p1_*.nc",
        "tos_day_bcc-csm1-1-m_historical_r1i1p1_*.nc",
        "tos_day_BNU-ESM_historical_r1i1p1_*.nc",
        "tos_day_CanESM2_historical_r1i1p1_*.nc",
        "tos_day_CCSM4_historical_r6i1p1_*.nc",
        "tos_day_CNRM-CM5_historical_r1i1p1_*.nc",
        "tos_day_FGOALS-g2_historical_r1i1p1_*.nc",
        "tos_day_GFDL-CM3_historical_r1i1p1_*.nc",
        "tos_day_GFDL-ESM2M_historical_r1i1p1_*.nc",
        "tos_day_IPSL-CM5A-MR_historical_r1i1p1_18500101-20051231.nc",
        "tos_day_GISS-E2-H_historical_r6i1p3_*.nc",
        "tos_day_MIROC5_historical_r1i1p1_*.nc",
        "tos_day_MIROC-ESM_historical_r1i1p1_*.nc",
        "tos_day_CSIRO-Mk3-6-0_historical_r1i1p1_*.nc",
        "tos_day_MRI-CGCM3_historical_r1i1p1_*.nc",
        "tos_day_NorESM1-M_historical_r1i1p1_*.nc"]
        
    for f1,f2 in zip(toslist,subsetlist):
        try:
            main(f1=f1,f2=f2,output_dir="subset/",inputvar="tos",outputvar="sst")
        except Exception as e:
            print(e,f1)