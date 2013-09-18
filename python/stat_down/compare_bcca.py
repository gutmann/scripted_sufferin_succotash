#!/usr/bin/env python
from glob import glob
import sys
import gc

import numpy as np
import matplotlib.pyplot as plt

import swim_io
import date_fun
from bunch import Bunch

bcca_search=["ccsm/*.nc"]
obs_search=["obs/*.nc"]
sar_search=["SAR/*.nc"]

bcca_search=["ccsm/*1950*.nc"]
obs_search=["obs/data.19[5-6]*.nc"]
sar_search=["SAR/BCSAR_Prec_6km_19[5-6]*.nc"]

bcca_search=["ccsm/*1970*.nc","ccsm/*1990*.nc"]
obs_search=["obs/data.19[7-9]*.nc","obs/data.200[0-5]*.nc"]
sar_search=["SAR/BCSAR_Prec_6km_19[7-9]*.nc","SAR/BCSAR_Prec_6km_200[0-5]*.nc"]


def wet_frac(data):
    """Calculate wet day fraction over the first array axis"""
    # size of one time slice
    nwet=np.zeros(data.shape[1:])
    # threshold to use for wet day fraction
    threshold=0
    # loop over all times
    for i in range(data.shape[0]):
        # find wet days in the current time period
        tmp=np.where(data[i,...]>threshold)
        # if there were more than 0 wet days
        if len(tmp[0])>0:
            # add one to all grid cells with a wet day
            nwet[tmp]+=1
    
    # calculate wet day fraction as the number of wet days / total number of days
    return nwet/float(data.shape[0])

def calc_stats(data,mask=None):
    """calculate basic statistics, annual means, Wet-day fraction, (99%,other...)
    
    should use stat_down.stats
    in most (all?) datasets, missing values are 1e20, in case they are 9999 in the future use a lower threshold to mask the data
    """
    # data structure to store output in
    stats=Bunch()
    
    # if we weren't supplied a mask calculate it from the current data
    if mask==None:
        mask=data[0,...]>9000
    # store mean annual precip in a masked array
    stats.mean=np.ma.array(np.mean(data,axis=0)*365.25,mask=mask)
    print(data.shape)
    # stddev takes a while to calculate and we don't care that much, add it later if desired
    # stdval=np.ma.array(np.std(data,axis=0),mask=mask)
    
    # store wet day fraction as a masked array
    stats.wetfrac=np.ma.array(wet_frac(data),mask=mask)

    return stats
    
def load_stats(filesearch,varname):
    """load varname from a host of files subset to a reasonable region"""
    # find the filenames to workwith
    filenames=glob(filesearch[0])
    # if there was more than one search pattern specified, search for those patterns too
    for i in range(1,len(filesearch)):
        filenames.extend(glob(filesearch[i]))
    # glob doesn't sort filenames correctly, we just want them sorted alphabetically (which = chronologically)
    filenames.sort()
    
    # read the first file in
    cur_data=swim_io.read_nc(filenames[0],varname).data
    
    # assuming monthly data (most variability in data length) calculate a total time length that should be an overestimate
    time_length=len(filenames)*31
    # now look at the length of data in the current/first file and if it was annual (or decadal...) assume that value for future files instead
    if len(filenames)*cur_data.shape[0]>time_length:
        # +10 is just an arbitrary number so we are less likely to over run again (assuming e.g. 10 more leap years per file)
        time_length=len(filenames)*(cur_data.shape[0]+10)
    # pre-allocate data arrays to minimize memory creation/destruction
    outputdata=np.zeros((time_length,cur_data.shape[1],314),dtype=np.float32)
    # current "end of data" pointer
    last=cur_data.shape[0]
    # store current data into the output data array
    outputdata[:last,...]=cur_data[:,:,:314]
    # this is a counter to keep track of files remaining in case we need to calculate a new output array
    remainingfiles=len(filenames)-1
    # loop over remaining files
    for f in filenames[1:]:
        cur_data=swim_io.read_nc(filenames[0],varname).data
        curlength=cur_data.shape[0]
        # if the data will overflow our output array increase the size of the output array
        if (last+curlength)>outputdata.shape[0]:
            print("    allocating more memory...")
            sys.stdout.flush()
            # create a new output array with a larger size
            # +10 is just an arbitrary number so we are less likely to over run again (assuming e.g. 10 more leap years per file)
            newoutputdata=np.zeros((last+(curlength+10)*remainingfiles,outputdata.shape[1],outputdata.shape[2]))
            # copy over old data
            newoutputdata[:last,...]=outputdata
            # delete old outputdata
            del outputdata
            # point outputdata to the new array
            outputdata=newoutputdata
            # make sure garbage collection runs
            gc.collect()
        # keeps track of how many files we have left to process to make it easy to calculate a reasonable new output array size if necessary
        remainingfiles-=1
        # add current data to the output data
        outputdata[last:last+curlength,...]=cur_data[:,:,:314]
        # update the end of data "pointer"
        last+=curlength
    
    # calculate statistics for this data set (e.g. mean, stdev, 99%, wetfrac,...)
    print("    Calculating Stats")
    sys.stdout.flush()
    stats=calc_stats(outputdata[:last,...])
    swim_io.write(filenames[0].split("/")[-1][:5]+"_mean",stats.mean)
    swim_io.write(filenames[0].split("/")[-1][:5]+"_wetfrac",stats.wetfrac)
    # once we have stats, make sure we clean up the memory from this dataset
    del outputdata
    gc.collect()
    
    return stats
        
        
    
def vis_stats(stats,filename="annual.png"):
    """visualize statistics"""
    plt.figure(figsize=(10,10))
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
    sys.stdout.flush()
    bcca_stats=load_stats(bcca_search,"pr")
    print("Loading Obs")
    sys.stdout.flush()
    obs=load_stats(obs_search,"Prec")
    print("Loading AR-BCCA")
    sys.stdout.flush()
    sar=load_stats(sar_search,"Prec")
    
    # not valid anymore, need to save from within load_stats when we want to deal with e.g. monthly data
    # dates=date_fun.mjd2date(np.arange(bcca.shape[0])+date_fun.date2mjd(1940,1,1,0,0))
    
    # store stats in a single data structure for the visualization routine
    stats=Bunch(bcca=bcca_stats,obs=obs,sar=sar)
    print("Visualizing...")
    sys.stdout.flush()
    vis_stats(stats)

if __name__ == '__main__':
    main()