#!/usr/bin/env python
from __future__ import print_function
import glob
import os,sys,re
import traceback

import numpy as np
from scipy import stats
plt=None
# import matplotlib.pyplot as plt

import mygis

def read_geosubset(filename,geo):
    geodata = mygis.read_geo(filename)
    lat,lon = geodata.lat, geodata.lon
    goodlat = np.where((lat>=geo[0]) & (lat<=geo[1]))[0]
    goodlon = np.where((lon>=geo[2]) & (lon<=geo[3]))[0]
    subset  = [goodlat[0],goodlat[-1]+1,goodlon[0],goodlon[-1]+1]
    return subset


def load_data(filesearch="*pr.19*.nc",variable="pr",n=None):
    files=glob.glob(filesearch)
    files.sort()
    geosubset=[34,43.9,360-115,360-100]
    geosubset=[25,52.7,360-124.7,360-67] #CONUS
    subset=read_geosubset(files[0],geosubset)
    # print(files[0],subset)
    # load in ALL the data (note this is likely to be 4-17GB)
    # on data vis clusters (128GB RAM) this is fine, 
    # on personal machines (4-8GB RAM) it's tricky or a a non-starter
    d=[]
    print("Loading "+str(len(files))+" files.")
    d=mygis.read_files(filesearch,variable,axis=0)
    d=d[:,subset[0]:subset[1],subset[2]:subset[3]]
    return d

if plt!=None:
    def vis_data(rs,minlag,maxlag,timelags,outputfile=""):
        for lag in range(minlag,maxlag+1):
            plt.clf();
            plt.imshow(rs[lag-1,:,:],vmin=0.3,vmax=1,origin="lower")
            plt.title("Spatial Lag="+str(lag))
            plt.colorbar()
            plt.draw()
            plt.savefig(outputfile+"Spatial_"+str(lag)+".png")

        for lag in range(1,timelags):
            plt.clf();
            plt.imshow(rs[lag-1+maxlag,:,:],vmin=0,vmax=0.7,origin="lower")
            plt.title("Temporal Lag="+str(lag))
            plt.colorbar()
            plt.draw()
            plt.savefig(outputfile+"Temporal_"+str(lag)+".png",dpi=150)


def main(d=None,outputfile="auto_correlations"):
    if d==None:
        d=load_data()
    shape=d.shape
    minlag=1
    maxlag=50
    timelags=4
    rs=np.zeros((maxlag-minlag+1+timelags+1,shape[1],shape[2]))+1E20
    delta=(shape[1]-maxlag*2)/20.0
    current=0.0
    for i in range(shape[1]):
        if (i-maxlag)>=current:
            current+=delta
            print(str(int(np.round(100.0*(i)/(shape[1]-1))))+"%",end=" ")
            sys.stdout.flush()
        for j in range(shape[2]):
            if d[0,i,j]<1e10:
                for lag in range(minlag,min(shape[1]-i-1,shape[2]-j-1,maxlag+1)):
                    r=0.0
                    n=0.0
                    # check for non-fill locations
                    if (d[0,i+lag,j]<1E10):
                        r2,p=stats.pearsonr(d[:,i,j],d[:,i+lag,j])
                        r+=r2
                        n+=1
                    if d[0,i,j+lag]<1E10:
                        r4,p=stats.pearsonr(d[:,i,j],d[:,i,j+lag])
                        r+=r4
                        n+=1
                    
                    if n>0:
                        rs[lag-1,i,j]=r/n
                    
                for t in range(1,timelags):
                    r,p=stats.pearsonr(d[t:,i,j],d[:-t,i,j])
                    rs[maxlag+t,i,j]=r
    print("...got data.")
    if glob.glob(outputfile+".nc"):
        os.rename(outputfile+".nc",outputfile+"_old.nc")
    mygis.write(outputfile,rs)
    
    if plt!=None:
        vis_data(rs,minlag,maxlag,timelags, outputfile)
    

def obs_driver():
    datasets=["../obs/maurer.125","../obs/uw.0625"]
    variables=["pr","tasmax","tasmin"]
    for d in datasets:
        for v in variables:
            try:
                filesearch="/".join([d, v, "*"+v+"*200[1-8].nc"])
                data=load_data(filesearch,v)
                outputname="autocorrelation_"+d.split("/")[-1]+"_"+v
                main(data,outputname)
                
            except Exception,e:
                print("While working on "+d+" "+v)
                print("Unexpected exception"+str(e))
                traceback.print_exc()                

def driver():
    # methods=["SDmon","CAe","SDe","CAe0","CAe1","SDe0","SDe1","CA","SARe0","SARe1"]
             # "CAcold","CAdry","SDcold","SDdry","CAhot","CAwet","SDhot","SDwet"]
    BCs=["BC",""]
    models=["ncep","narr"]
    resolutions=["12km","6km"]
    variables=["pr","tasmax","tasmin"]
    methods=["CA","SDmon_c","SD","SAR"]
    
    for b in BCs:
        for forc in models:
            for r in resolutions:
                for v in variables:
                    for m in methods:
                        print(m,forc,r,v,b)
                        extra=""
                        if forc=="ncep":extra="*gauss" #differentiate NCEP T2m from sig995 for temperature
                        
                        filesearch="/".join([m, forc, v, b+m[:2]+"*"+r+extra+"*200[1-8]*.nc"])
                        files=glob.glob(filesearch)
                        # print(filesearch)
                        if len(files)>5:
                            try:
                                d=load_data(filesearch,v,n=None)
                                output="-".join([m,forc,v,b+r])
                                main(d,output+"autocorrelation")
                                
                            except Exception,e:
                                print("While working on :"+"-".join([m,forc,v,b+r]))
                                print("Unexpected exception"+str(e))
                                traceback.print_exc()
                            

if __name__ == '__main__':
    # driver()
    obs_driver()
    # main()
