#!/usr/bin/env python
import os,glob,copy
import fnmatch
import numpy as np
import sys

import mygis
from bunch import Bunch
from bilin_regrid import regrid,load_geoLUT, vinterp

geographic_vars = ["lat","lon","latitude", "longitude", "XLAT","XLONG", "XLAT_M","XLONG_M","time"] #,"z","PH","PHB"]

def main(f1=None, f2=None, zvar=None, zbasevar=None,skip_z0=False,
            zrefvar=None, zrefbasevar=None, timevar=None,
            output_dir="./", inputvar=None, outputvar=None,
            regrid_vertical=True, regrid_horizontal=True):
    """docstring for main"""
    if f1!=None:
        inputfile=f1
    else:
        inputfile="tos_day_ACCESS1-3_historical_r1i1p1_19500101-19591231.nc"
    if f2!=None:
        reference=f2
    else:
        reference="../../completed/historical/access/subset/ta_6hrLev_ACCESS1-3_historical_r1i1p1_1950010106-1951010100.nc"

    # find all files if more than one was specified.
    files=glob.glob(inputfile)
    files.sort()

    if regrid_vertical and regrid_horizontal:
        if zvar!=None:
            geographic_vars.append(zvar)
            if zbasevar!=None:
                geographic_vars.append(zbasevar)

    if timevar!=None:
        geographic_vars.append(timevar)

    print("Reading geo data")
    geo1 = mygis.read_geo(files[0])
    geo2 = mygis.read_geo(reference)
    while len(geo1.lat.shape)>2:
        geo1.lat=geo1.lat[0]
        geo1.lon=geo1.lon[0]

    while len(geo2.lat.shape)>2:
        geo2.lat=geo2.lat[0]
        geo2.lon=geo2.lon[0]

    if geo1.latatts is None:
        latatts=Bunch(standard_name="latitude",long_name="latitude coordinate",units="degrees_north")
    else:
        latatts=geo1.latatts
    if geo1.lonatts is None:
        lonatts=Bunch(standard_name="longitude",long_name="longitude coordinate",units="degrees_east")
    else:
        lonatts=geo1.lonatts

    master_vars=[Bunch(data=geo2.lat,name="lat",dims=("lat","lon"),dtype="f",attributes=latatts),
                 Bunch(data=geo2.lon,name="lon",dims=("lat","lon"),dtype="f",attributes=lonatts)]
    if timevar:
        timeatts=Bunch(standard_name="time")
        master_vars.append(Bunch(data=None,name=timevar,dims=("time",),dtype="d",attributes=timeatts))

    missing_value=0
    # find the variable names that we need to regrid
    # assume we are regridding everything that is not a "geographic" variable (e.g. lat, lon, time)
    if inputvar is None:
        fdata = mygis.Dataset(files[0])
        varlist = []
        for v in fdata.variables:
            if (not (str(v) in geographic_vars)):
                varlist.append(str(v))
        fdata.close()
    else:
        varlist=[inputvar]

    if outputvar is None:
        outputvar=varlist

    print(varlist)
    print(latatts)
    print(lonatts)
    if timevar:print(timeatts)
    # ncdata=mygis.read_nc(files[0],inputvar,returnNCvar=True)
    if regrid_horizontal:
        geoLUT = load_geoLUT(lat1=geo1.lat,lon1=geo1.lon,lat2=geo2.lat,lon2=geo2.lon, winhalfsize=5)#,
                        # mask=(ncdata.data[0,...]==missing_value),winhalfsize=7)
        print(geoLUT.shape)
    # ncdata.ncfile.close()
    if zrefvar:
        zout=mygis.read_nc(reference,zrefvar).data
        if zrefbasevar!=None:
            zout+=mygis.read_nc(reference,zrefbasevar).data

    for f in files:
        regridded_z=None
        if zvar:
            zin=mygis.read_nc(f,zvar).data
            if zbasevar!=None:
                zin+=mygis.read_nc(f,zbasevar).data
            # if zvar=="GHT":
            #     zin/=9.81
            if skip_z0:
                zin=zin[:,1:]

        extra_vars = copy.deepcopy(master_vars)
        global_atts= mygis.read_atts(f,global_atts=True)

        if timevar:
            time_atts=mygis.read_atts(f,varname=timevar)
            extra_vars[2].data=mygis.read_nc(f,timevar).data
            extra_vars[2].attributes=time_atts

        for i,v in enumerate(varlist):
            print("Reading data: "+f+" for variable:"+v)
            data=mygis.read_nc(f,v).data
            # the first model layer in erai is on the surface, but the second layer is at 1000mb (often below the surface)
            #  this lets you skip the surface layer and just interpolate over the pressure levels.
            if skip_z0:
                if len(data.shape)==4:
                    data=data[:,1:]

            data_atts=mygis.read_atts(f,varname=v)
            if len(data.shape)==4:
                dims=("time","lev","lat","lon")
            elif len(data.shape)==3:
                dims=("time","lat","lon")
            elif len(data.shape)==2:
                dims=("lat","lon")
            else:
                print("Error with variable "+v)
                print("can't geointerpolate 1D variable, skipping")
                dims=mygis.read_dims(f,varname=v)
            print(mygis.read_dims(f,varname=v))
            print(dims)

            print("Regridding")
            if (regrid_horizontal and (len(data.shape)>1)):
                print(data.shape)
                output=regrid(data,geoLUT=geoLUT)
            else:
                output=data

            if zvar:
                if regrid_vertical:
                    if regridded_z is None:
                        if regrid_horizontal:
                            regridded_z = regrid(zin,geoLUT=geoLUT)
                        else:
                            regridded_z = zin
                    if len(output.shape)==4:
                        nz = zout.shape[1]
                        for t in range(output.shape[0]):
                            output[t,:nz] = vinterp(output[t],inputz=regridded_z[t],outputz=zout[t*4+2])
                        output=output[:,:nz,:,:]
            print(data.shape)
            print(output.shape)

            if v!=varlist[-1]:
                extra_vars.append(Bunch(data=output,name=v,dims=dims,dtype="f",attributes=data_atts))
            else:
                print("Writing data")
                outputfile=f.split("/")[-1]
                print(outputfile)
                print(output_dir)
                print(v)
                print(outputvar[i])
                outputfile=output_dir+"regrid_"+outputfile.replace(v,outputvar[i])
                print(outputfile)
                print(output.shape)
                for e in extra_vars:
                    print(e.name, e.dims, e.data.shape)
                if outputfile[2:]!=inputfile:
                    mygis.write(outputfile,output,dims=dims,varname=outputvar[i],attributes=data_atts,
                                global_attributes=global_atts,extravars=extra_vars)#, format="NETCDF4")
                else:
                    print("Error, outputfile==inputfile")

if __name__ == '__main__':

    # need to set up and argparse interface
    if len(sys.argv)==3:
        main(f1=sys.argv[1],f2=sys.argv[2])#,output_dir=sys.argv[3],inputvar=sys.argv[4],outputvar=sys.argv[4])
    else:
        main(f1=sys.argv[1],f2=sys.argv[2],output_dir=sys.argv[3],inputvar=sys.argv[4],outputvar=sys.argv[5])

    # various one offs for specific datasets
    # if not os.path.exists("subset"):
    #     os.mkdir("subset")

    # convert (rotated pole?) grid sea surface temperatures (tos) to rectilinear grid (ta)
    # subsetlist=[
    #         "../../completed/historical/access/subset/ta_6hrLev_ACCESS1-3_historical_r1i1p1_1950010106-1951010100.nc",
    #         "../../completed/historical/bcc/subset/ta_6hrLev_bcc-csm1-1-m_historical_r1i1p1_195001010000-195003141800.nc",
    #         "../../completed/historical/bnu/subset/ta_6hrLev_BNU-ESM_historical_r1i1p1_1950010100-2005123118.nc",
    #         "../../completed/historical/canesm/subset/ta_6hrLev_CanESM2_historical_r1i1p1_195001010000-195012311800.nc",
    #         "../../completed/historical/ccsm/subset/ta_6hrLev_CCSM4_historical_r6i1p1_1950010106-1950033118.nc",
    #         "../../completed/historical/cnrm_cm5/subset/ta_6hrLev_CNRM-CM5_historical_r1i1p1_195001010600-195002010000.nc",
    #         "../../completed/historical/fgoals/subset/ta_6hrLev_FGOALS-g2_historical_r1i1p1_1950010106-1951010100.nc",
    #         "../../completed/historical/gfdl_cm3/subset/ta_6hrLev_GFDL-CM3_historical_r1i1p1_1950010100-1950123123.nc",
    #         "../../completed/historical/gfdl_esm/subset/ta_6hrLev_GFDL-ESM2M_historical_r1i1p1_1951010100-1955123123.nc",
    #         "../../completed/historical/ipsl_mr/subset/ta_6hrLev_IPSL-CM5A-MR_historical_r1i1p1_1950010103-1959123121.nc",
    #         "../../completed/historical/giss_e2h/subset/ta_6hrLev_GISS-E2-H_historical_r6i1p3_195001010600-195007010000.nc",
    #         "../../completed/historical/miroc5/subset/ta_6hrLev_MIROC5_historical_r1i1p1_1950010100-1950013118.nc",
    #         "../../completed/historical/miroc_esm/subset/ta_6hrLev_MIROC-ESM_historical_r1i1p1_1950010106-1950020100.nc",
    #         "../../completed/historical/mk3/subset/ta_6hrLev_CSIRO-Mk3-6-0_historical_r1i1p1_195001010600-195101010000.nc",
    #         "../../completed/historical/mri_cgcm3/subset/ta_6hrLev_MRI-CGCM3_historical_r1i1p1_195001010000-195001311800.nc",
    #         "../../completed/historical/noresm/subset/ta_6hrLev_NorESM1-M_historical_r1i1p1_1950010100-1950063018.nc"]
    #
    # toslist=[
    #     "tos_day_ACCESS1-3_historical_r1i1p1_*.nc",
    #     "tos_day_bcc-csm1-1-m_historical_r1i1p1_*.nc",
    #     "tos_day_BNU-ESM_historical_r1i1p1_*.nc",
    #     "tos_day_CanESM2_historical_r1i1p1_*.nc",
    #     "tos_day_CCSM4_historical_r6i1p1_*.nc",
    #     "tos_day_CNRM-CM5_historical_r1i1p1_*.nc",
    #     "tos_day_FGOALS-g2_historical_r1i1p1_*.nc",
    #     "tos_day_GFDL-CM3_historical_r1i1p1_*.nc",
    #     "tos_day_GFDL-ESM2M_historical_r1i1p1_*.nc",
    #     "tos_day_IPSL-CM5A-MR_historical_r1i1p1_18500101-20051231.nc",
    #     "tos_day_GISS-E2-H_historical_r6i1p3_*.nc",
    #     "tos_day_MIROC5_historical_r1i1p1_*.nc",
    #     "tos_day_MIROC-ESM_historical_r1i1p1_*.nc",
    #     "tos_day_CSIRO-Mk3-6-0_historical_r1i1p1_*.nc",
    #     "tos_day_MRI-CGCM3_historical_r1i1p1_*.nc",
    #     "tos_day_NorESM1-M_historical_r1i1p1_*.nc"]
    #
    # for f1,f2 in zip(toslist,subsetlist):
    #     try:
    #         main(f1=f1,f2=f2,output_dir="subset/",inputvar="tos",outputvar="sst")
    #     except Exception as e:
    #         print(e,f1)
    # subsetlist=glob.glob("pr_*_ncep_1979*.nc")
    # subsetlist.sort()
    # newlist=[]
    # for i in range(len(subsetlist)):
    #     if fnmatch.fnmatch(subsetlist[i],"pr_MM5I_ncep_*.nc"):
    #         pass
    #     else:
    #         newlist.append(subsetlist[i])
    #
    # for f1 in newlist:
    #     print(f1)
    #     try:
    #         main(f1=f1,f2="pr_MM5I_ncep_1991010103.nc",output_dir="subset/",inputvar="pr",outputvar="pr")
    #     except Exception as e:
    #         print(e,f1)
    #
