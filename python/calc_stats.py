from __future__ import print_function
import argparse
import sys

import numpy as np
from scipy import stats

import swim_io
import stats_driver
    
def wetdry(data,threshold=0):
    """Calculate wet/dry statistics for data given a precip threshold
    
    data = 3D grid of precip data [time x lat x lon]
    threshold = precip threshold to use as definition of a wet day
    returns wet day fraction, mean wet spell length, and mean dry spell length
        on a 2D grid [lat x lon]
    """
    wetdays=np.zeros(data.shape)
    wetdays[data>threshold]=1
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
    
    wetspell_length=total_wetdays/wetspells
    dryspell_length=(data.shape[0]-total_wetdays)/dryspells
    
    return(wetfraction,wetspell_length,dryspell_length)
    
def mean(data,nyears=9):
    """Calculate a mean annual total from a daily time series spanning nyears"""
    return(data.sum(axis=0)/nyears)

def interannual(data,yearstarts):
    """Calculate interannual variability in mean annual values
    
    data = 3D grid [time x lat x lon]
    yearstarts = list of indices into the time dimension specifying 
                the starting point for each year
    """
    annshape=list(data.shape)
    annshape[0]=len(yearstarts)
    annual_vals=np.zeros(annshape)
    for i in range(annshape[0]):
        startpt=yearstarts[i]
        if i==annshape[0]-1:
            endpt=None
        else:
            endpt=yearstarts[i+1]
        annual_vals[i,...]=np.mean(data[startpt:endpt,...],axis=0)
    
    return(np.std(annual_vals,axis=0))

def histogram(data,precip=True):
    """Calculate a histogram for data
    
    For precip calculate on a log scale. 
    returns histgram bin value and center point as a tuple
    """
    if precip:
        usedata=np.log10(data[data>0])
        minmax=[-1,3]
    else:
        usedata=data
        minmax=[-30,60]
    
    v,x=np.histogram(data,bins=30,range=minmax)
    
    if precip:
        x=10.0**x
    
    # return bin centers instead of edges
    x=(x[1:]+x[:-1])/2
    return(v,x)
    
def calc_fritch(data):
    """Calculate the fritch index from tmin data (C)"""
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
    """Calculate various temperature indices from tmin data
    
    data = 3D grid of tmin (C or K) [time, lat, lon]
    yearstarts = list of indices into the time dimension specifying year starting points
    returns mean annual Frost days per year and Fritch index
    on a 2D grid [lat, lon]
    """
    if data.min()>100:
        data-=273.15
        datawerekelvin=True
    
    # n frost days in time series
    frost_days=np.choose(data<=0,(0,1))
    # frost days per year = nfrost days / nyears
    frost_dpy=frost_days.sum(axis=0)/float(len(yearstarts))
    
    fritchshape=list(data.shape)
    fritchshape[0]=len(yearstarts)
    fritch_indicies=np.zeros(fritchshape)

    # loop over years calculating the fritch index for each year
    for i in range(fritchshape[0]):
        startpt=yearstarts[i]
        if i==fritchshape[0]-1:
            endpt=None
        else:
            endpt=yearstarts[i+1]
        fritch_indicies[i,...]=calc_fritch(data[startpt:endpt,...])
        
    return(frost_dpy,fritch_indicies.mean(axis=0))
    
def p99(data):
    """Calculate the 99th and 1st percentile for data along the first index
    
    data can be a 2D or 3D array, 
        sorts data on axis 0 and 
        returns the 99th and 1st percentile on a 1D or 2D grid
    """
    sorted_data=np.sort(data,axis=0)
    return (sorted_data[np.round(data.shape[0]*0.99),...],
            sorted_data[np.round(data.shape[0]*0.01),...])


def calc_extreme_value(params,distribution,nyear,datafraction=1.0):
    """Given parameters for and a distribution, 
    calculate the n-year return interval storm
    datafraction is the fraction of input data that were 
        used when calculating the distribution"""
    probability=1.0/(nyear*365.25/datafraction)
    stepsize=10
    x=np.arange(-1000,1000,stepsize).astype('d')
    pdf=distribution.pdf(x,*params)
    pos=np.where((pdf>=probability) & (pdf>0))[0]
    
    if len(pos)==0:return -9999
    pos=pos[-1]

    x=np.arange(x[pos]-2*stepsize,x[pos]+2*stepsize,0.1).astype('d')
    pdf=distribution.pdf(x,*params)
    pos=np.where((pdf<=probability) &(pdf>0))[0]
    
    # if x[pos[0]]<0:return 0
    return x[pos[0]]


def extremes(data, dist_name="gamma",verbose=True):
    """Calculate the 2, 5, 10 and 50yr return intervals for data using a specified distribution
    
    data = 3D grid of e.g. precip data [time, lat, lon]
    dist_name = the distribution to use when fitting data (gamma, weibull, or exponential)
    """
    # set numpy error handling functions to ignore invalid data because we may get some, 
    # but store the old settings so we can restore them after this function
    old_settings=np.seterr(invalid='ignore')
    
    # select the distribution to use
    distribution=stats.gamma
    if dist_name=="weibull":
        distribution=stats.weibull_min
    elif dist_name=="exponential":
        distribution=stats.expon

    # return intervals to calculate storms for
    year_intervals=[2,10,50,100]
    shape=data.shape
    extremes=np.zeros((len(year_intervals),shape[1],shape[2]))
    if verbose:print("Calculating Extremes")
    # loop over lat dimension
    for i in range(shape[1]):
        # this may take a while so output progress to the screen if verbose output was specified
        if verbose:
            print("progress= {0}/{1}      \r".format(i,shape[1]),end="")
            sys.stdout.flush()
        # loop over lon dimension
        for j in range(shape[2]):
            curdata=data[:,i,j]
            # just use the top half of the non-zero values in data
            medianval=np.median(curdata[curdata>0])
            usevals=np.where(curdata>=medianval)[0]
            if len(usevals)>0:
                # assuming there were some data to fit, fit the data
                params=distribution.fit(curdata[usevals])
                # now generate a list of calculate extremes for each return interval
                curextremes=[calc_extreme_value(params,distribution,
                                                year,len(usevals)/float(shape[0])) 
                                                for year in year_intervals]
                # and store the data in an output array
                extremes[:,i,j]=np.array(curextremes)
    print("\nFinished")
    # restore error checking stats
    new_settings=np.seterr(invalid=old_settings["invalid"])
    return extremes
    
    
    