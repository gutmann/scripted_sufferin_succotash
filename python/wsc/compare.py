import datetime
from glob import glob

import numpy as np
from bunch import Bunch

import regrid, snodas, recon, wrf


comparison_domain=Bunch(lat=[35.5,42],lon=[-109,-104.5])

def load_wrf(years):
    all_data=wrf.load("wrf/SWE_daily.nc")
    wrfdata=[]
    for y in years:
        curstart=np.where(all_data.dates>=datetime.datetime(y,1,1,0,0))[0][0]
        curend=np.where(all_data.dates>=datetime.datetime(y,12,31,0,0))[0][0]
        wrfdata.append(Bunch(data=all_data.data[curstart:curend+1,:,:],
                            lat=all_data.lat,lon=all_data.lon,dates=all_data.dates[curstart:curend+1]))
    return wrfdata

def all():
    """Compare all SWE products"""
    years=range(2004,2009)
    
    print("Loading WRF data")
    wrfdata=load_wrf(years)
    junk=wrf.stats(wrfdata[0].data)
    snodasdata=[]
    recondata=[]
    print("Loading recon and snodas data")
    for y in years:
        print("   "+str(y))
        snodasfile="snodas/SWE_Daily0600UTC_WesternUS_"+str(y)+".dat"
        snodasdata.append(snodas.load(snodasfile,startyear=y))
        
        reconfile="SWE_SNODIS/"+str(y)+"/swe.dat"
        recondata.append(recon.load(reconfile,startyear=y))


    geo=wrfdata[0]
    latsub=np.where((geo.lat>=comparison_domain.lat[0]) & (geo.lat<=comparison_domain.lat[1]))[0]
    latsub=[latsub.min(),latsub.max()]
    lonsub=np.where((geo.lon>=comparison_domain.lon[0]) & (geo.lon<=comparison_domain.lon[1]))[1]
    lonsub=[lonsub.min(),lonsub.max()]

    print(comparison_domain)
    print(latsub)
    print(lonsub)
    
    r_lut=None
    s_lut=None
    rstat=[]
    sstat=[]
    wstat=[]
    print("computing Stats")
    for s,r,w in zip(snodasdata,recondata,wrfdata):
        print(w.dates[0])
        r_lut,rdata=regrid.agg(r,w.lat[latsub[0]:latsub[1],lonsub[0]:lonsub[1]],w.lon[latsub[0]:latsub[1],lonsub[0]:lonsub[1]],geo_lut=r_lut)
        s_lut,sdata=regrid.agg(s,w.lat[latsub[0]:latsub[1],lonsub[0]:lonsub[1]],w.lon[latsub[0]:latsub[1],lonsub[0]:lonsub[1]],geo_lut=s_lut)
        
        rstat.append(recon.stats(rdata.data))
        sstat.append(snodas.stats(sdata.data))
        wstat.append(wrf.stats(w.data[:,latsub[0]:latsub[1],lonsub[0]:lonsub[1]]))
    
    
    return (rstat,sstat,wstat)
