#!/usr/bin/env python
import glob,sys,os

import numpy as np

import mygis
import bilin_regrid
from bunch import Bunch

gcm_varnames=["hus","ps","rh","ta","ua","va","z","sst","slp"]
erai_varnames=["Q_GDS4_ISBL_S123","LNSP_GDS4_HYBL_S123","R_GDS4_ISBL_S123","T_GDS4_ISBL_S123",
                "U_GDS4_ISBL_S123","V_GDS4_ISBL_S123","Z_GDS4_ISBL_S123","SSTK_GDS4_SFC_S123","MSL_GDS4_SFC_S123"]
erai_fileprefixes=["SC","PS","SC","SC","UV","UV","SC","SF","SF"]

erai_name="{prefix}_month{month:02}_mean.nc"
global_outputfile="regridded_ERAi_to_{gcm}_month{month:02}.nc"
gcm_filename="{gcm}_mean_month{month:02}_{varname}.nc"

dim3d=("p","lat","lon")
dim2d=("lat","lon")
dims=[dim3d]*9
dims[1]=dim2d
dims[7]=dim2d
dims[8]=dim2d
latdim=dim2d
londim=dim2d
pdim=dim3d
eraid="erai/" #directory containing ERAi monthly mean data

atts=[Bunch(long_name="specific humidity",units="kg/kg"),
      Bunch(long_name="surface pressure",units="Pa"),
      Bunch(long_name="relative humidity",units="%"),
      Bunch(long_name="air temperature",units="K"),
      Bunch(long_name="eastward wind",units="m/s"),
      Bunch(long_name="northward wind",units="m/s"),
      Bunch(long_name="geopotential height",units="m"),
      Bunch(long_name="sea surface temperature",units="K"),
      Bunch(long_name="sea level pressure",units="Pa")
      ]

patts=Bunch(long_name="pressure",units="Pa")
latatts=Bunch(long_name="latitude",units="degrees north")
lonatts=Bunch(long_name="longitude",units="degrees east")

def read_gcm_geo(filename,pfile):
    """docstring for read_geo"""
    geo=mygis.read_geo(filename)
    lat=geo.lat
    lon=geo.lon
    if lon.max()>180:
        lon=lon-360
    if os.path.isfile(pfile):
        p=mygis.read_nc(pfile,"pressure").data
    else:
        p=None
    return Bunch(p=p,lat=lat,lon=lon)

def read_erai_geo(filename):
    """docstring for read_geo"""
    geo=mygis.read_geo(filename)
    lat=geo.lat
    lon=geo.lon
    if lon.max()>180:
        lon=lon-360
    p=mygis.read_nc(filename,"p").data
    return Bunch(p=p,lat=lat,lon=lon)


def write_output(outputfile,data,geo):
    """docstring for write_output"""
    vardata=[
        Bunch(data=data[erai_varnames[i]],name=v,dtype="f",attributes=atts[i],dims=dims[i])
        for i,v in enumerate(gcm_varnames)
    ]
    # for v in vardata:
    #     print(v.name,v.dims,v.data.shape)
    vardata.insert(0,Bunch(data=geo.lon,name="lon",dtype="f",attributes=lonatts,dims=londim))
    vardata.insert(0,Bunch(data=geo.lat,name="lat",dtype="f",attributes=latatts,dims=latdim))
    # if geo.p!=None:
    #     vardata.append(Bunch(data=geo.p,name="p",dtype="f",attributes=patts,dims=pdim))
    
    if len(geo.p.shape)==1:
        ny,nx=geo.lat.shape
        geo.p=geo.p[:,np.newaxis,np.newaxis].repeat(ny,axis=1).repeat(nx,axis=2)
    
    mygis.write(outputfile,data=geo.p,varname="p",dtype="f",attributes=patts,dims=pdim,
                history="Regridding performed by cmip/era2gcm.py",extravars=vardata)

def main(gcm="ccsm"):
    """docstring for main"""
    gcmd=gcm+"/" #directory containing GCM monthly mean data
    erai_geo=read_erai_geo(eraid+"SC_month01_mean.nc")
    gcm_geo=read_gcm_geo(gcmd+"geo.nc",gcmd+gcm+"_month01_mean_plevel.nc")
    print("Creating Geographic Look up table")
    geoLUT=bilin_regrid.load_geoLUT(erai_geo.lat,erai_geo.lon,gcm_geo.lat,gcm_geo.lon)
    for month in range(1,13):
        print("Month: "+str(month))
        gcm_geo=read_gcm_geo(gcmd+"geo.nc",gcmd+gcm+"_month{0:02}_mean_plevel.nc".format(month))
        if gcm_geo.p!=None:
            print("Creating Vertical Look up table")
            pLUT=bilin_regrid.load_vLUT(erai_geo.p*100,gcm_geo.p)
        
        output_era_data=dict()
        for i,(p,v) in enumerate(zip(erai_fileprefixes,erai_varnames)):
            print("loading:"+v)
            raw_era_data=mygis.read_nc(eraid+erai_name.format(prefix=p,month=month),v).data
            print("Interpolating")
            if (v[:4]=="LNSP") or (v[:3]=="MSL") or (v[:4]=="SSTK"):
                output_era_data[v]=bilin_regrid.regrid(raw_era_data[np.newaxis,:,:],geoLUT=geoLUT)[0,...]
            else:
                era_on_gcm_grid=bilin_regrid.regrid(raw_era_data,geoLUT=geoLUT)
                if gcm_geo.p!=None:
                    output_era_data[v]=bilin_regrid.vinterp(era_on_gcm_grid,vLUT=pLUT)
                else:
                    output_era_data[v]=era_on_gcm_grid
            if v=="Z_GDS4_ISBL_S123":
                output_era_data[v]/=9.8
            # gcm_data=mygis.read_nc(gcmd+gcm_filename.format(gcm=gcm,month=month,varname=gcm_varnames[i]),gcm_varnames[i]).data
        if gcm_geo.p==None:
            gcm_geo.p=erai_geo.p*100
        print("Writing data")
        write_output(global_outputfile.format(month=month,gcm=gcm),output_era_data,gcm_geo)
        
    
        

if __name__ == '__main__':
    gcm="ccsm"
    if len(sys.argv)>1:
        gcm=sys.argv[1]
    main(gcm=gcm)