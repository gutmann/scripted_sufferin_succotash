from __future__ import print_function
import sys
import numpy as np
from bunch import Bunch

def load_geoLUT(lat1,lon1,lat2,lon2):
    N1=lat1.shape
    N2=lat2.shape
    geoLUT=np.empty((N2[0],N2[1],2),dtype=list)
    for i in range(N1[0]):
        print(i,N1[0])
        for j in range(N1[1]):
            dists=(lat1[i,j]-lat2)**2+(lon1[i,j]-lon2)**2
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            if not geoLUT[newy,newx,0]:
                geoLUT[newy,newx,0]=list()
                geoLUT[newy,newx,1]=list()
            geoLUT[newy,newx,0].append(i)
            geoLUT[newy,newx,1].append(j)
    return geoLUT

def regrid_hi2low(data,lat1=None,lon1=None,lat2=None,lon2=None,geoLUT=None,FillValue=1E20):
    """docstring for regrid_hi2low"""
    if geoLUT==None:
        if len(lat2.shape)==1:
            lon2,lat2=np.meshgrid(lon2,lat2)
        if len(lat1.shape)==1:
            lon1,lat1=np.meshgrid(lon1,lat1)
        print("Computing LUT")
        geoLUT=load_geoLUT(lat1,lon1,lat2,lon2)
    
    twoD=False
    if len(data.shape)==2:
        twoD=True
        data=data[np.newaxis,:,:]
    N=data.shape
    outputdata=np.empty((N[0],geoLUT.shape[0],geoLUT.shape[1]))
    N2=outputdata.shape
    print("processing...")
    # NOTE, this might be speed up substantially by inlining C code... 
    # except that geoLUT is an array of variable length lists... 
    # not sure how weave handles that would probably have to convert to a larger array first
    for i in range(N2[0]):
        if (i%25)==0:
            print(i,end=" ")
            sys.stdout.flush()
        for j in range(N2[1]):
            for k in range(N2[2]):
                if geoLUT[j,k,0]:
                    curdata=data[i,geoLUT[j,k,0],geoLUT[j,k,1]]
                    tmp=np.where(curdata<1E15)
                    if len(tmp[0])>0:
                        outputdata[i,j,k]=curdata[curdata<1E15].mean()
                    else:
                        outputdata[i,j,k]=FillValue
                        
                else:
                    outputdata[i,j,k]=FillValue
                
    if twoD:
        outputdata=outputdata.reshape(outputdata.shape[1:])
            
    return Bunch(data=outputdata,geo=geoLUT)