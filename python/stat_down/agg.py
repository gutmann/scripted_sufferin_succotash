import numpy as np
from bunch import Bunch
import mygis as swim_io
import sys


def dataXn(data,nsamples):
    """aggregate an NxM data array into an nsamples x nsamples array
    
    output=dataXn(data,nsamples)
    
    NOTE: if nsamples<0 it is treated as a stride length instead"""
    ny,nx=data.shape
    if nsamples>0:
        ystep=np.floor(ny/float(nsamples))
        xstep=np.floor(nx/float(nsamples))
        xsamples=nsamples
        ysamples=nsamples
    elif nsamples<0:
        xstep=-nsamples
        ystep=-nsamples
        xsamples=np.floor(nx/float(xstep))
        ysamples=np.floor(ny/float(ystep))
        
    outputdata=np.zeros((ysamples,xsamples))
    ystart=0
    for i in range(ysamples):
        xstart=0
        for j in range(xsamples):
            curdata=data[ystart:ystart+ystep,xstart:xstart+xstep]
            outputdata[i,j]=np.mean(curdata[np.isfinite(curdata)])
            xstart+=xstep
        ystart+=ystep
    return outputdata
    
def data2huc(data,hucmask,huc,minvalue=-100,maxvalue=1E5):
    """Compute mean of data where hucmask==huc
    
        if data is a 3D array loop over the first index and return a time series
        if data is a 2D array return the mean value
        
        also eliminates values <min or >max values (default=-100,1E5)
        
        usage: output=data2huc(data,hucmask,huc,minvalue=-100,maxvalue=1E5)
    """
    if len(data.shape)>2:
        # if data is a 3D array, loop over the first index... not ideal
        outputdata=np.zeros(data.shape[0])
        for i in range(data.shape[0]):
            curdata=data[i,...]
            curdata=curdata[hucmask==huc]
            gooddata=np.where((curdata>=minvalue) & (curdata<1E5) & np.isfinite(curdata))
            if len(gooddata[0])>0:
                outputdata[i]=np.mean(curdata[gooddata])
            else:
                outputdata[i]=-9999
        return outputdata
    # if data is just a 2D array, simply return a mean value
    else:
        curdata=data[hucmask==huc]
        return np.mean(curdata[(curdata>=0) & (curdata<1E5) & np.isfinite(curdata)])

def data2hucs(data,hucmask,minvalue=-100,maxvalue=1E5):
    """Compute the mean of data for each huc in hucmask
    
    output=data2hucs(data,hucmask,minvalue=0,maxvalue=1E5)
        data = np.array(dims= 2 or 3)
        hucmask=np.array(dims=2) matching last 2 dims of data
        min/max values = valid range of values in data
        
    returns    a list of either mean values (if data==2D)
            or a list of np.arrays of len data.shape[0] (if data==3D)
    """
    huclist=np.unique(hucmask)
    huclist=huclist[huclist>0]
    hucdata=[]
    
    for huc in huclist:
        hucdata.append(data2huc(data,hucmask,huc,minvalue=minvalue,maxvalue=maxvalue))
        
    return Bunch(data=hucdata,hucs=huclist)
    
