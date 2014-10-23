#!/usr/bin/env python
import sys
import numpy as np
import mygis

from bunch import Bunch
import glob

filesearch_dict=dict(uv="ei.moda.an.pl/ncfiles/ei.moda.an.pl.regn128uv.????{0:02}*.nc",
                     sc="ei.moda.an.pl/ncfiles/ei.moda.an.pl.regn128sc.????{0:02}*.nc",
                     ml="ei.moda.an.ml/ncfiles/ei.moda.an.ml.regn128sc.????{0:02}*.nc")

pl_suf="_GDS4_ISBL_S123"
ml_suf="_GDS4_HYBL_S123"
varlist_dict=dict(uv=["U"+pl_suf,"V"+pl_suf],
                  sc=["Z"+pl_suf,"Q"+pl_suf,"R"+pl_suf,"T"+pl_suf],
                  ml=["LNSP"+ml_suf])
             
geovars_dict=dict(uv=Bunch(lat="g4_lat_1",lon="g4_lon_2",p="lv_ISBL0"),
                  sc=Bunch(lat="g4_lat_1",lon="g4_lon_2",p="lv_ISBL0"),
                  ml=Bunch(lat="g4_lat_0",lon="g4_lon_1"))

outputfile_dict=dict(uv="UV_month{0:02}_mean.nc",
                     sc="SC_month{0:02}_mean.nc",
                     ml="PS_month{0:02}_mean.nc")

global_atts=dict()
global_atts["U"+pl_suf]=dict(long_name="Eastward Wind",units="m s**-1")
global_atts["V"+pl_suf]=dict(long_name="Northward Wind",units="m s**-1")
global_atts["Z"+pl_suf]=dict(long_name="Geopotential Height",units="m**2 s**-2")
global_atts["Q"+pl_suf]=dict(long_name="Specific Humidity",units="kg kg**-1")
global_atts["R"+pl_suf]=dict(long_name="Relative Humidity",units="%")
global_atts["T"+pl_suf]=dict(long_name="Temperature",units="K")
global_atts["LNSP"+ml_suf]=dict(long_name="Surface Pressure",units="hPa")
global_atts["lat"]=dict(long_name="latitude",units="degrees_north")
global_atts["lon"]=dict(long_name="longitude",units="degrees_east")
global_atts["p"]=dict(long_name="pressure",units="hPa")

dim3d=("p","lat","lon")
dim2d=("lat","lon")
latdim=("lat",)
londim=("lon",)
dim1d=("p",)
global_dims=dict()
global_dims["U"+pl_suf]=dim3d
global_dims["V"+pl_suf]=dim3d
global_dims["Z"+pl_suf]=dim3d
global_dims["Q"+pl_suf]=dim3d
global_dims["R"+pl_suf]=dim3d
global_dims["T"+pl_suf]=dim3d
global_dims["LNSP"+ml_suf]=dim2d
global_dims["lat"]=latdim
global_dims["lon"]=londim
global_dims["p"]=dim1d


def write_outputfile(filename,data,lat,lon,p=None):
    """write the output file"""
    mainvar=data.keys()[0]
    extravars=[]
    
    main_atts=global_atts[mainvar]
    main_dims=global_dims[mainvar]
    
    for varname in data.keys()[1:]:
        extravars.append(Bunch(data=data[varname],name=varname,dtype="f",
                              attributes=global_atts[varname],dims=global_dims[varname]))
    
    extravars.append(Bunch(data=lat,name="lat",dtype="f",attributes=global_atts["lat"],dims=global_dims["lat"]))
    extravars.append(Bunch(data=lon,name="lon",dtype="f",attributes=global_atts["lon"],dims=global_dims["lon"]))
    if p!=None:
        extravars.append(Bunch(data=p,name="p",dtype="i",attributes=global_atts["p"],dims=global_dims["p"]))
    
    mygis.write(filename,data[mainvar],varname=mainvar,dtype="f",
                dims=main_dims,attributes=main_atts,extravars=extravars)


def main(etype="ml"):
    """docstring for main"""
    
    geovars=geovars_dict[etype]
    varlist=varlist_dict[etype]
    outputfile=outputfile_dict[etype]
    filesearch=filesearch_dict[etype]

    for month in range(1,13):
        print("Processing month: {}".format(month))
        
        #find all the files to process for this month
        files=glob.glob(filesearch.format(month))
        files.sort() # not strictly necessary
        
        # dict to store outputdata in
        outputdata=dict()
        #loop over all variables
        for v in varlist:
            print("   "+v)
            n=0 #count the number of files we have loaded
            # loop over all files
            for f in files:
                #read in the current variable for the current file
                curdata=mygis.read_nc(f,v).data
                #if it is LNSP (Natural log of surface pressure) convert to hPa
                if v[:4]=="LNSP":
                    curdata=np.exp(curdata)
                
                # if these data exist in the outputdata dictionary add the current data
                # and add one to the count (n)
                if v in outputdata:
                    outputdata[v]+=curdata
                    n+=1
                # if these data do not exist in the outputdata dictionary store the current data
                # and set the count (n) to 1
                else:
                    outputdata[v]=curdata
                    n=1
                    
            # divide by n to get the mean for the current variable
            # n SHOULD equal the number of files
            if n!=len(files):
                print("Warning not all files had variable:"+v)
            if n!=0:
                outputdata[v]/=float(n)
            else:
                print("No files with this variable:"+v)
                print("Example file:"+files[i])
        
        # load geographic data:
        lat=mygis.read_nc(f,geovars.lat).data
        lon=mygis.read_nc(f,geovars.lon).data
        # try to load pressure level data (valid for all but the LNSP model level data)
        try:
            p=mygis.read_nc(f,geovars.p).data
        except:
            p=None
        
        #write the output data for the current month
        write_outputfile(outputfile.format(month),outputdata,lat=lat,lon=lon,p=p)
    

if __name__ == '__main__':
    print("""Set var to run:
             sc  ml  uv""")
    if len(sys.argv)>1:
        etype=sys.argv[1]
    else:
        etype="ml"
    main(etype)