#!/usr/bin/env python
from glob import glob
import gc

import numpy as np
import matplotlib.pyplot as plt

import swim_io
import date_fun
from bunch import Bunch

bcca_search=["ccsm/*.nc"]
obs_search=["obs/*.nc"]
sar_search=["SAR/*.nc"]

def wet_frac(data):
    """Calculate wet day fraction over the first array axis"""
    nwet=np.zeros(data.shape[1:])
    for i in range(data.shape[0]):
        tmp=np.where(data[i,...]>threshold)
        if len(tmp[0])>0:
            nwet[tmp]+=1
    
    return nwet/float(data.shape[0])

def calc_stats(data,mask=None):
    """calculate basic statistics, annual means, Wet-day fraction, (99%,other...)
    
    should use stat_down.stats
    """
    stats=Bunch()
    
    if mask==None:
        mask=data[0,...]>9000
    stats.mean=np.ma.array(np.mean(data,axis=0)*365.25,mask=mask)
    # stdval=np.ma.array(np.std(data,axis=0),mask=mask)
    stats.wetfrac=np.ma.array(wet_frac(data,axis=0),mask=mask)

    return stats
    
def load_stats(filesearch,varname):
    """load varname from a host of files subset to a reasonable region"""
    filenames=glob(filesearch[0])
    for i in range(1,len(filesearch)):
        filenames.extend(filesearch[i])
    filenames.sort()
    d1=swim_io.read_nc(filenames[0],varname).data
    outputdata=np.zeros((d1.shape[0],d1.shape[1],314),dtype=np.float32)
    outputdata[:]=d1[:,:,:314]
    for f in filenames[1:]:
        d1=swim_io.read_nc(filenames[0],varname).data
        newoutputdata=np.zeros((outputdata.shape[0]+d1.shape[0],outputdata.shape[1],outputdata.shape[2]),dtype=np.float32)
        newoutputdata[:outputdata.shape[0],...]=outputdata
        newoutputdata[outputdata.shape[0]:,...]=d1[:,:,:314]
        outputdata=newoutputdata
        gc.collect()
    
    stats=calc_stats(outputdata)
    del outputdata
    gc.collect()
    return stats
        
        
    
def vis_stats(stats,filename="annual.png"):
    """visualize statistics"""
    plt.figure(xs=6,ys=6)
    plt.subplot(231)
    plt.imshow(stats.obs.mean,vmax=3000)
    plt.subplot(232)
    plt.imshow(stats.sar.mean,vmax=3000)
    plt.subplot(233)
    plt.imshow(stats.bcca.mean,vmax=3000)

    plt.subplot(234)
    plt.imshow(stats.obs.mean-stats.obs.mean,vmax=500,vmin=-500)
    plt.subplot(235)
    plt.imshow(stats.sar.mean-stats.obs.mean,vmax=500,vmin=-500)
    plt.subplot(236)
    plt.imshow(stats.bcca.mean-stats.obs.mean,vmax=500,vmin=-500)
    
    plt.savefig("means_"+filename)

    plt.clf()
    plt.subplot(231)
    plt.imshow(stats.obs.wetfrac,vmax=1)
    plt.subplot(232)
    plt.imshow(stats.sar.wetfrac,vmax=1)
    plt.subplot(233)
    plt.imshow(stats.bcca.wetfrac,vmax=1)
    plt.subplot(234)
    plt.imshow(stats.obs.wetfrac-stats.obs.wetfrac,vmax=1,vmin=-1)
    plt.subplot(235)
    plt.imshow(stats.sar.wetfrac-stats.obs.wetfrac,vmax=1,vmin=-1)
    plt.subplot(236)
    plt.imshow(stats.bcca.wetfrac-stats.obs.wetfrac,vmax=1,vmin=-1)
    
    plt.savefig("wetfrac_"+filename)
    

def main():
    """"Compare BCCA, ARBBCA, obs"""
    print("Loading BCCA")
    bcca_stats=load_stats(bcca_search,"pr")
    print("Loading Obs")
    obs=load_stats(obs_search,"Prec")
    print("Loading AR-BCCA")
    sar=load_stats(sar_search,"Prec")
    
    # dates=date_fun.mjd2date(np.arange(bcca.shape[0])+date_fun.date2mjd(1940,1,1,0,0))
    
    stats=Bunch(bcca=bcca_stats,obs=obs,sar=sar)
    # print("Calculating Statistics")
    # stats=calc_stats(bcca,obs,sar)
    print("Visualizing...")
    vis_stats(stats)

if __name__ == '__main__':
    main()