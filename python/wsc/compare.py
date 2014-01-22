#!/usr/bin/env python
import datetime
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from bunch import Bunch

import regrid, snodas, recon, wrf, modscag


comparison_domain=Bunch(lat=[35.5,42],lon=[-109,-104.5])

def split2years(all_data,years):
    wrfdata=[]
    for y in years:
        curstart=np.where(all_data.dates>=datetime.datetime(y,1,1,0,0))[0][0]
        curend=np.where(all_data.dates>=datetime.datetime(y,12,31,0,0))[0][0]
        wrfdata.append(Bunch(data=all_data.data[curstart:curend+1,:,:],
                            lat=all_data.lat,lon=all_data.lon,dates=all_data.dates[curstart:curend+1]))
    return wrfdata

def load_wrf(years):
    all_data=wrf.load("wrf/SWE_daily.nc")
    return split2years(all_data,years)

def sca():
    """Compare SCA between MODSCAG and WRF
    
    Assumes a directory structure of 
        wrf/SWE_daily.nc
        MODSCAG/fsca2008.dat
                dateselect.txt
    """
    print("Loading WRF data")
    wrfdata=load_wrf([2008])
    print("Loading MODSCAG data")
    modsca=modscag.load("MODSCAG/fsca2008.dat")
    
    geo=wrfdata[0]
    latsub=np.where((geo.lat>=comparison_domain.lat[0]) & (geo.lat<=comparison_domain.lat[1]))[0]
    latsub=[latsub.min(),latsub.max()]
    lonsub=np.where((geo.lon>=comparison_domain.lon[0]) & (geo.lon<=comparison_domain.lon[1]))[1]
    lonsub=[lonsub.min(),lonsub.max()]
    
    print("Regridding...")
    sca_lut=None
    sca_lut,scadata=regrid.agg(modsca,
                                geo.lat[latsub[0]:latsub[1],lonsub[0]:lonsub[1]],
                                geo.lon[latsub[0]:latsub[1],lonsub[0]:lonsub[1]],
                                geo_lut=sca_lut)
                                
    wsca=wrf.sca(wrfdata[0].data[:,latsub[0]:latsub[1],lonsub[0]:lonsub[1]]/1000.0)
    
    print("visualizing good images")
    for i in range(len(modsca.gooddates.indices)):
        plt.clf()
        plt.subplot(121)
        plt.imshow(wsca[modsca.gooddates.indices[i],:,:],vmax=1.1)
        plt.title(modsca.gooddates.dates[i]);
        plt.subplot(122)
        plt.imshow(scadata.data[modsca.gooddates.indices[i],:,:],vmax=1.1)
        plt.title(modsca.gooddates.dates[i]);
        plt.draw()
        plt.savefig("SCA_comparison_{0:03}.png".format(i))

    print("visualizing all images")
    for i in range(wsca.shape[0]):
        plt.clf()
        plt.subplot(121)
        plt.imshow(wsca[i,:,:],vmax=1.1)
        plt.title(modsca.dates[i]);
        plt.subplot(122)
        plt.imshow(scadata.data[i,:,:],vmax=1.1)
        plt.title(modsca.dates[i]);
        plt.draw()
        plt.savefig("complete_SCA_comparison_{0:03}.png".format(i))

    
    z=swim_io.read_nc("wrf/4km_wrf_output.nc","HGT").data
    zsub=z[0,latsub[0]:latsub[1],lonsub[0]:lonsub[1]]

    print("visualizing sca=f(z)")
    thresholds=[0.1,0.5,0.9]
    for t1 in thresholds:
        wrf_mid_z=np.zeros(len(modsca.gooddates.indices))
        sca_mid_z=np.zeros(len(modsca.gooddates.indices))
        for i in range(len(modsca.gooddates.indices)):
            wrf_pts=np.where((wsca[modsca.gooddates.indices[i],:,:]>=t1-0.05)
                    &(wsca[modsca.gooddates.indices[i],:,:]<=t1+0.05))
            if len(wrf_pts[0])>0:
                wrf_mid_z[i]=zsub[wrf_pts].mean()
            sca_pts=np.where((sca_grid[modsca.gooddates.indices[i],:,:]>=t1-0.05)
                    &(sca_grid[modsca.gooddates.indices[i],:,:]<=t1+0.05))
            if len(sca_pts[0])>0:
                sca_mid_z[i]=zsub[sca_pts].mean()
    
        plt.clf()
        plt.plot(modsca.gooddates.dates,wrf_mid_z,'.',label="WRF")
        plt.plot(modsca.gooddates.dates,sca_mid_z,'.',label="MODSCAG")
        plt.title("Mean elevation at {0}% SCA".format(t1*100))
        plt.legend(loc=4)
        plt.ylim(1500,4000)
        plt.draw()
        plt.savefig("z_sca_{}.png".format(t1*100))
    
        
            
    
    

def all():
    """Compare all SWE products"""
    years=range(2004,2009)
    
    print("Loading WRF data")
    wrfdata=load_wrf(years)
    # junk=wrf.stats(wrfdata[0].data)
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

if __name__ == '__main__':
    sca()

# for i in range(len(years)):
#     plt.clf();
#     subplot(1,2,1);plt.imshow(sstat[i].melt_date[:,17:]-rstat[i].melt_date[:,:-17],cmap=plt.cm.seismic,vmin=-50,vmax=50)
#     plt.title("SNODAS - Noahmp")
#     plt.colorbar();
#     plt.title("Melt Date - SNODAS - Noahmp")
#     subplot(1,2,2);plt.imshow(sstat[0].melt_date[:,17:]-wstat[0].melt_date[:-11,17:-1],cmap=plt.cm.seismic,vmin=-50,vmax=50)
#     plt.title("Melt Date - SNODAS - Noah")
#     plt.colorbar();
#     plt.savefig("melt-"+str(years[i])+".png")
#     
# for i in range(len(years)):
#     plt.clf();
#     subplot(1,2,1);plt.imshow(sstat[0].peak[:,17:]*1000-rstat[0].peak[:,:-17],cmap=plt.cm.seismic,vmin=-600,vmax=600)
#     plt.title("Peak SWE - SNODAS - Noahmp")
#     plt.colorbar();
#     subplot(1,2,2);plt.imshow(sstat[0].peak[:,17:]*1000-wstat[0].peak[:-11,17:-1],cmap=plt.cm.seismic,vmin=-600,vmax=600)
#     plt.title("Peak SWE - SNODAS - Noah")
#     plt.colorbar();
#     plt.savefig("peak-"+str(years[i])+".png")
#     
