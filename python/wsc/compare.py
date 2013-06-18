from glob import glob
import numpy as np

from bunch import Bunch

import regrid, snodas, recon, wrf


comparison_domain=Bunch(lat=[35.5,42],lon=[-109,-104.5])

def all():
    """Compare all SWE products"""
    years=range(2004,2009)
    snodasdata=[]
    recondata=[]
    for y in years:
        snodasfile="snodas/SWE_Daily0600UTC_WesternUS_"+str(y)+".dat"
        snodasdata.append(snodas.load(snodasfile,startyear=y))
        
        reconfile="SWE_SNODIS/"+str(y)+"/swe.dat"
        recondata.append(recon.load(reconfile,startyear=y))

    geo=snodasdata[0]
    latsub=np.where((geo.lat>=comparison_domain.lat[0]) & (geo.lat<=comparison_domain.lat[1]))[0]
    latsub=[latsub.min(),latsub.max()]
    lonsub=np.where((geo.lon>=comparison_domain.lon[0]) & (geo.lon<=comparison_domain.lon[1]))[0]
    lonsub=[lonsub.min(),lonsub.max()]
    
    rlut=None
    rstat=[]
    sstat=[]
    for s,r in zip(snodasdata,recondata):
        rlut,thisdata=regrid.agg(r,s.lat[latsub[0]:latsub[1]],s.lon[lonsub[0]:lonsub[1]],geo_lut=rlut)
        
        rstat.append(recon.stats(thisdata.data))
        sstat.append(snodas.stats(s.data[:,latsub[0]:latsub[1],lonsub[0]:lonsub[1]]))
        
    return (rstat,sstat)
