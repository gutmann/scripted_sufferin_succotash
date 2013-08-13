#!/usr/bin/env python
from glob import glob

import numpy as np
import matplotlib.pyplot as plt

import swim_io
import date_fun
from bunch import Bunch

bcca_search="ccsm/*.nc"
obs_search="obs/*.nc"
sar_search="SAR/*.nc"

def wet_frac(data):
    """Calculate wet day fraction over the first array axis"""
    nwet=np.zeros(data.shape[1:])
    for i in range(data.shape[0]):
        tmp=np.where(data[i,...]>threshold)
        if len(tmp[0])>0:
            nwet[tmp]+=1
    
    return nwet/float(data.shape[0])

def calc_stats(bcca,obs,sar):
    """calculate basic statistics comparing to obs in annual means, Wet-day fraction, 99%"""
    bcca_stats=Bunch(); obs_stats=Bunch(); sar_stats=Bunch()
    
    mask=obs[0,...]>9000
    bcca_stats.mean=np.ma.array(np.mean(bcca,axis=0)*365.25,mask=mask)
    obs_stats.mean=np.ma.array(np.mean(obs,axis=0)*365.25,mask=mask)
    sar_stats.mean=np.ma.array(np.mean(sar,axis=0)*365.25,mask=mask)

    # bcca_stats.std=np.ma.array(np.std(bcca,axis=0),mask=mask)
    # obs_stats.std=np.ma.array(np.std(obs,axis=0),mask=mask)
    # sar_stats.std=np.ma.array(np.std(sar,axis=0),mask=mask)
    
    bcca_stats.wetfrac=np.ma.array(wet_frac(bcca),mask=mask)
    obs_stats.wetfrac=np.ma.array(wet_frac(obs),mask=mask)
    sar_stats.wetfrac=np.ma.array(wet_frac(sar),mask=mask)
    
    return Bunch(bcca=bcca_stats,obs=obs_stats,sar=sar_stats)
    
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
    bcca=np.concatenate(swim_io.read_files(bcca_search,"pr"))
    print("Loading Obs")
    obs=np.concatenate(swim_io.read_files(obs_search,"Prec"))
    print("Loading AR-BCCA")
    sar=np.concatenate(swim_io.read_files(sar_search,"Prec"))
    dates=date_fun.mjd2date(np.arange(bcca.shape[0])+date_fun.date2mjd(1940,1,1,0,0))
    
    print("Calculating Statistics")
    stats=calc_stats(bcca,obs,sar)
    print("Visualizing...")
    vis_stats(stats)

if __name__ == '__main__':
    main()