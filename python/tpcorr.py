#!/usr/bin/env python
import traceback,os
import numpy as np
from scipy import stats
import re

import matplotlib.pyplot as plt
from matplotlib import cm

import swim_io
# import stats_driver
from stat_down import driver as stats_driver
from bunch import Bunch

def find_var(filename,possibleNames):
    """loads data from a file when the variable name is one of a list of possible names
    
    This could move into swim_io.py...?
    """
    outputvar=None
    try:
        while outputvar==None:
            curname=possibleNames.pop()
            try:
                data=swim_io.read_nc(filename,curname).data
                return data
            except KeyError:
                pass
    except StopIteration:
        return None
    

def read_geo(filename):
    """reads geographic data (lat,lon) from a file
    
    This could move into swim_io.py...?
    """
    lat=find_var(filename,["lat","XLAT","latitude"])
    lon=find_var(filename,["lon","XLONG","longitude"])
    if lon.max()>180:
        lon=lon-360
    # if len(lat.shape)==1:
    #     lat,lon=np.meshgrid(lat,lon)
    return lat,lon

def getsubset(filename,geodata):
    """Find subset indices in a file given geographic restrictions in geodata"""
    lat,lon=read_geo(filename)
    latstart=np.where(lat>=geodata.latmin)[0][0]
    latend=np.where(lat<=geodata.latmax)[0][-1]
    lonstart=np.where(lon>=geodata.lonmin)[0][0]
    
    lonend=np.where(lon<=geodata.lonmax)[0][-1]
    if lonstart>lonend:
        lonend=np.where(lon>=geodata.lonmin)[0][-1]
        lonstart=np.where(lon<=geodata.lonmax)[0][0]
    return Bunch(latstart=latstart,latend=latend,lonstart=lonstart,lonend=lonend)

def load_all_data(files,variable,n=None,geosubset=None,skipaug03=True):
    """load in ALL the data (note this is likely to be 4-17GB)
    
    on data vis clusters (128GB RAM) this is fine, 
    on desktop machines (4-8GB RAM) it's tricky or a a non-starter
    """
    d=[]
    if geosubset:
        subset=getsubset(files[0],geosubset)
        print(subset)
    if variable=="pr":
        threshold=200
        fill=0
    else:
        curdata=swim_io.read_nc(files[0],variable,returnNCvar=True)
        if np.median(curdata.data[0:10,...])<200:
            # Temperature is in K
            threshold=350
            fill=290
        else:
            # Temperature is in C
            threshold=70
            fill=20
        curdata.ncfile.close()
            
    for f in files[:n]:
        if geosubset:
            filedata=swim_io.read_nc(f,variable,returnNCvar=True)
            curdata=filedata.data[:,subset.latstart:subset.latend+1,
                                    subset.lonstart:subset.lonend+1]
            tmp=np.where(curdata>threshold)
            if len(tmp[0])>0:
                curdata[tmp]=fill
            
            filedata.ncfile.close()
        else:
            curdata=swim_io.read_nc(f,variable).data
        # if skipaug03:
        #     if re.match(".*2003_08.*", f):
        #         curdata=None
        #     elif re.match(".*2003.*", f):
        #         print(f)
        #         if curdata.shape[0]>243:
        #             # remove august 2003, mainly a problem for the 4 and 6km datasets. 
        #             curdata=np.vstack([curdata[:211,...],curdata[243:,...]])
        if curdata!=None:
            d.append(curdata)
    d=np.concatenate(d,axis=0)
    return d


def write_output(filename,data):
    swim_io.write("pr_corr_"+filename,data[0])
    swim_io.write("pr_pvals_"+filename,data[1])
    
    img=np.ma.array(data[0])
    img.mask=np.array(img.shape,dtype=np.bool)
    img.mask[:]=True
    img.mask[data[1]<0.05]=False
    # img[:]=data[0]
    
    plt.clf()
    plt.imshow(img,origin="lower",vmin=-0.3,vmax=0.3,cmap=cm.seismic)
    plt.colorbar()
    plt.title(filename)
    plt.draw()
    plt.savefig("pr_corr_"+filename+".png")


def calc_correlation(d1,d2,n=None):
    sz=d1.shape
    rout=np.zeros(sz[1:])
    pout=np.zeros(sz[1:])
    if n==None:
        n=min(d1.shape[0],d2.shape[0])
    for i in range(sz[1]):
        for j in range(sz[2]):
            r1,p1=stats.pearsonr(d1[:n,i,j],d2[:n,i,j])
            rout[i,j]=r1
            pout[i,j]=p1
    return [rout,pout]

def remove_seasonal_cycle(data):
    """Removes the seasonal cycle from a long time series of data
    
    First averages across all years, then applies a +/- 30day 
    averaging window (that wraps around Jan. 1)
    """
    sz=list(data.shape)
    sz[0]=365
    annualmeans=np.zeros(sz)
    count=0
    # first average all years together
    for i in np.arange(0.0,data.shape[0],365.25):
        if (int(i)+365)<=data.shape[0]:
            annualmeans+=data[int(i):int(i)+365,...]
            count+=1
    annualmeans/=count
    windowsize=30
    smoothedmeans=np.zeros(sz)
    # then apply a +/-30 day smoothing window
    for i in range(365):
        if i<windowsize:
            bottom=np.vstack([annualmeans[i-windowsize:,...],annualmeans[:i,...]])
        else:
            bottom=annualmeans[i-windowsize:i,...]
        if i+windowsize>sz[0]:
            top=np.vstack([annualmeans[:(i+windowsize)-sz[0],...],annualmeans[i:,...]])
        else:
            top=annualmeans[i:i+windowsize,...]
        smoothedmeans[i,...]=(top.mean(axis=0)+bottom.mean(axis=0))/2
    smmean=smoothedmeans.mean(axis=0)
    
    # Because that smoothing will decrease the magnitude of the summer max/winter min
    #  adjust it by increasing the seasonal variability by 5% (rough hack that looked right for some data)
    adjusted=smoothedmeans-smmean
    adjusted*=1.05
    adjusted+=smmean
    
    # next remove that smoothed, adjusted annual cycle from the original data
    for i in np.arange(0.0,data.shape[0],365.25):
        if (int(i)+365)<=data.shape[0]:
            data[int(i):int(i)+365,...]-=adjusted
    return data

def tpcorr(infolist):
    # geosubset=Bunch(latmin=34.0,latmax=43.9,lonmin=245-360,lonmax=260-360)
    geosubset=Bunch(latmin=25.0,latmax=50.0,lonmin=-125.0,lonmax=-72.0)
    info=Bunch()
    for i in infolist:
        key=i.varname
        info[key]=i
    print("Loading data for:"+infolist[0].output)
    pr_data=load_all_data(info.pr.files,info.pr.varname,geosubset=geosubset)
    print("  tasmin...")
    tmin_data=load_all_data(info.tasmin.files,info.tasmin.varname,geosubset=geosubset)
    print("  tasmax...")
    tmax_data=load_all_data(info.tasmax.files,info.tasmax.varname,geosubset=geosubset)
    
    print("Removing temperature seasonal cycle")
    tmin_data=remove_seasonal_cycle(tmin_data)
    tmax_data=remove_seasonal_cycle(tmax_data)

    print("Computing correlations for:"+info.tasmax.output)
    ptmax_correlation=calc_correlation(pr_data,tmax_data)
    write_output(info.tasmax.output,ptmax_correlation)

    print("Computing correlations for:"+info.tasmin.output)
    ptmin_correlation=calc_correlation(pr_data,tmin_data)
    write_output(info.tasmin.output,ptmin_correlation)
    
def tpcorr_varcollector(files,varname,output,listing,extra):
    """Collects datasets until it has a complete set (tmin,tmax,pr) then processes
    
    Processing is performed by tpcorr (main function)
    This should probably be moved to stats_driver.py
    """
    # note: listing=[method,bc_status,forcing_model,resolution,variable]
    # print(extra)
    files.sort()
    if extra[0]==None:
        extra[0]=Bunch(files=files,varname=varname,output=output,listing=listing)
        if len(extra)==1:
            extra.extend([None,None])
    else:
        matchLastCall=True
        for this,last in zip(listing[:-1],extra[0].listing[:-1]):
            if this!=last:
                matchLastCall=False
        if matchLastCall:
            if extra[1]==None:
                extra[1]=Bunch(files=files,varname=varname,output=output,listing=listing)
            else:
                extra[2]=Bunch(files=files,varname=varname,output=output,listing=listing)
                tpcorr(extra)
                extra[0]=None
                extra[1]=None
                extra[2]=None
        else:
            print('Two successive calls did not match!:')
            print("   original="+str(extra[0].listing))
            print("   current ="+str(listing))
            extra[0]=Bunch(files=files,varname=varname,output=output,listing=listing)
            extra[1]=None
            extra[2]=None
    
if __name__ == '__main__':
    # collection=stats_driver.RunCollector()
    # stats_driver.drive(collection.collect)
    # tpcorr(collection)
    try:
        stats_driver.drive(tpcorr_varcollector,
                           yearsearch="200[0-5]",
                           resolutions=["12km"],
                           models=["ncep","narr"],
                           methods=["SDmon_c","SD","CA","SAR"],
                           # methods=["SDmon","CAe0","SDe0","CA","SARe0"],
                           BCs=["BC"],
                           stat=True,obs=True,runforce=True)
    except Exception as e:
        print("Error Unexpected Exception")
        print(e)
        traceback.print_exc()
        os._exit(1)
