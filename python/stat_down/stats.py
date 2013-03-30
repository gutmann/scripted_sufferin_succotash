#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys

import numpy as np
from scipy import stats

import swim_io
import stats_driver
    
def wetdry(data,threshold=0):
    wetdays=np.zeros(data.shape)
    for i in range(data.shape[0]):
        wetdays[i,data[i,...]>threshold]=1
    total_wetdays=wetdays.sum(axis=0)
    
    wetfraction=total_wetdays/data.shape[0]
    
    wetspells=np.zeros(data.shape[1:])
    dryspells=np.zeros(data.shape[1:])
    wetspells[wetdays[0,...]==1]=1
    dryspells[wetdays[0,...]==0]=1
    for i in range(1,data.shape[0]):
        new_wetdays=np.where((wetdays[i,...]==1) & (wetdays[i-1,...]==0))
        wetspells[new_wetdays]+=1

        # in the end n_dryspells~n_wetspells (+/- 1?)
        new_drydays=np.where((wetdays[i,...]==0) & (wetdays[i-1,...]==1))
        dryspells[new_drydays]+=1
    
    wetspells[wetspells==0]=1 #wetdays should be 0 so length will be 0/1
    wetspell_length=total_wetdays/wetspells
    dryspells[dryspells==0]=1 #drydays should be 0 so length will be 0/1
    dryspell_length=(data.shape[0]-total_wetdays)/dryspells
    
    return(wetfraction,wetspell_length,dryspell_length)
    
def mean(data,nyears=9,minval=-999,maxval=1e5):
    return(data.sum(axis=0)/nyears)

def interannual(data,yearstarts,fun=None):
    if fun==None:
        fun=np.mean
    annshape=list(data.shape)
    annshape[0]=len(yearstarts)
    annual_vals=np.zeros(annshape)
    for i in range(annshape[0]):
        startpt=yearstarts[i]
        if i==annshape[0]-1:
            endpt=None
        else:
            endpt=yearstarts[i+1]
        annual_vals[i,...]=fun(data[startpt:endpt,...],axis=0)
    
    return(np.std(annual_vals,axis=0))

def histogram(data,precip=True):
    if precip:
        usedata=[]
        for i in range(data.shape[0]):
            usedata.append(np.log10(data[i,data[i,...]>0]))
        usedata=np.concatenate(usedata)
        minmax=[-2,3]
    else:
        usedata=data
        minmax=[-30,60]
    
    v,x=np.histogram(usedata,bins=50,range=minmax)
    
    if precip:
        x=10.0**x
    
    # return bin centers instead of edges
    x=(x[1:]+x[:-1])/2
    return(v,x)
    
def calc_fritch(data):
    growing_threshold=5
    
    startdate=np.zeros(data.shape[1:])+999
    enddate=np.zeros(data.shape[1:])
    curwarmdays=np.zeros(data.shape[1:])
    
    for doy in range(data.shape[0]):
        warmdays=np.where(data[doy,...]>growing_threshold)
        colddays=np.where(data[doy,...]<=growing_threshold)
        
        if len(warmdays[0])>0:
            curwarmdays[warmdays]+=1
        if len(colddays[0])>0:
            curwarmdays[colddays]=0
        
        growing=np.where(curwarmdays==5)
        if len(growing[0])>0:
            startdate[growing]=np.choose((doy-5)<startdate[growing],(startdate[growing],doy-5))
            enddate[growing]=np.choose(doy>enddate[growing],(enddate[growing],doy))
        
    growing_season=enddate-startdate
    growing_season[growing_season<0]=0
    return growing_season
        
    
def temperature_indicies(data,yearstarts):
    if data.min()>100:
        data-=273.15
        datawerekelvin=True
    frost_days=np.choose(data<=0,(0,1))
    frost_dpy=frost_days.sum(axis=0)/float(len(yearstarts))
    
    fritchshape=list(data.shape)
    fritchshape[0]=len(yearstarts)
    fritch_indicies=np.zeros(fritchshape)

    for i in range(fritchshape[0]):
        startpt=yearstarts[i]
        if i==fritchshape[0]-1:
            endpt=None
        else:
            endpt=yearstarts[i+1]
        fritch_indicies[i,...]=calc_fritch(data[startpt:endpt,...])
        
    return(frost_dpy,fritch_indicies.mean(axis=0))
    
def p99(data):
    sorted_data=np.sort(data,axis=0)
    return (sorted_data[np.round(data.shape[0]*0.99),...],
            sorted_data[np.round(data.shape[0]*0.01),...])


def calc_extreme_value(params,distribution,nyear,datafraction=1.0):
    """Given a scipy distribution and associated parameters, 
    calculate the n-year return interval storm
    datafraction is the fraction of input data that were 
        used when calculating the distribution."""
    probability=1.0/(nyear*365.25/datafraction)
    return distribution.ppf(1-probability,*params)


def extremes(data, dist_name="gamma",verbose=True):
    if dist_name=="weibull":
        distribution=stats.weibull_min
    elif dist_name=="exponential":
        distribution=stats.expon
    elif dist_name=="gamma":
        distribution=stats.gamma
    else:
        return None
    old_settings=np.seterr(invalid='ignore')
        
    year_intervals=[2,10,50,100]
    shape=data.shape
    extremes=np.zeros((len(year_intervals),shape[1],shape[2]))
    if verbose:print("Calculating Extremes")
    for i in range(shape[1]):
        if verbose:
            print("progress= {0}/{1}      \r".format(i,shape[1]),end="")
            sys.stdout.flush()
        for j in range(shape[2]):
            curdata=data[:,i,j]
            medianval=np.median(curdata[curdata>0])
            # medianval=np.median(curdata[curdata>medianval])
            usevals=np.where(curdata>=medianval)[0]
            if len(usevals)>0:
                params=distribution.fit(curdata[usevals])
                curextremes=[calc_extreme_value(params,distribution,
                                                year,len(usevals)/float(shape[0])) 
                                                for year in year_intervals]
                extremes[:,i,j]=np.array(curextremes)
                if extremes[:,i,j].min()==0:print(i,j)
                
    print("\nFinished")
    new_settings=np.seterr(invalid=old_settings["invalid"])
    return extremes
    
    
    