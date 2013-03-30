#!/usr/bin/env python
from __future__ import print_function
import sys,os,traceback
import argparse

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

import swim_io
from nc_reader import NC_Reader
import stats_driver
import agg
from bunch import Bunch

def simple_load_data(files,varname):
    outputdata=[]
    for f in files:
        outputdata.append(swim_io.read_nc(f,varname).data)
    return np.concatenate(outputdata,axis=0)


def load_data(files,varname,extra):
    if extra[2]==None:
        return simple_load_data(files,varname)
    
    reader=NC_Reader(None,filelist=files,
                    geo_subset=extra[2],
                    readvars=[varname],ntimes=365,
                    latvar="lat",lonvar="lon")
    outputdata=[]
    nsamples=int(extra[3])
    hucfilename=extra[5]
    if hucfilename:
            hucreader=NC_Reader(None,filelist=[hucfilename],
                            geo_subset=extra[2],
                            readvars=["data"],ntimes=1,
                            latvar="lat",lonvar="lon")
            hucdata=hucreader.next()[0]
            hucreader.close()
        
    huclist=[]
    for data in reader:
        if nsamples==0:
            if hucfilename:
                aggdata=agg.data2hucs(data[0],hucdata,minvalue=0,maxvalue=500)
                data=np.array(aggdata.data).reshape((1,1,-1))
                huclist=aggdata.hucs
            outputdata.append(data)
        else:
            outputdata.append([agg.dataXn(data[0],nsamples)])
        
    reader.close()
    outputdata=np.concatenate(outputdata,axis=0)
    print(outputdata.shape)
    return Bunch(data=outputdata,hucs=huclist)

def calc_extreme_value(params,distribution,nyear,datafraction=1.0):
    """Given parameters for and a distribution, 
    calculate the n-year return interval storm
    datafraction is the fraction of input data that were 
        used when calculating the distribution"""
    probability=1.0/(nyear*365.0*datafraction)
    
    x=np.arange(0,10000,10)
    pdf=distribution.pdf(x,*params)
    pos=np.where((pdf>=probability) & (pdf>0))[0]
    
    if len(pos)==0:return 0
    pos=pos[-1]

    x=np.arange(x[pos]-10,x[pos]+10,0.1)
    pdf=distribution.pdf(x,*params)
    pos=np.where((pdf<=probability) &(pdf>0))[0]

    if x[pos[0]]<0:return 0
    return x[pos[0]]
    

def calc_extremes(files,v,output_base,info,extra,data=None,huclist=[]):
    print("Loading Data")
    sys.stdout.flush()
    if data==None:
        data=load_data(files,v,extra)
        huclist=data.hucs
        data=data.data
    nt,ny,nx=data.shape
    
    print(output_base)
    
    distribution=extra[0]
    outputdir=extra[4]
    outputname=outputdir+"extremes_"+output_base
    maxvalue=1000.0
    params=[]
    extreme_vals=[]
    
    nday_sum=1
    year_intervals=[2,10,50,100]
    nday_intervals=range(1,6)
    curparams=[0,0,0,0]
    curextremes=range(len(year_intervals))
    shortcut=1
    print("Calculating extreme values")
    for i in range(0,ny,shortcut):
        print("progress = "+str(float(i)/ny*100)+"%      ",end="\r")
        sys.stdout.flush()
        params.append([])
        extreme_vals.append([])
        for j in range(0,nx,shortcut):
            threshold=extra[1]*np.median(data[data[:,i,j]<maxvalue,i,j])
            usevals=np.where((data[:,i,j]>threshold)&(data[:,i,j]<(maxvalue*nday_sum)))[0]
            if len(usevals)>0:
                curparams=distribution.fit(data[usevals,i,j])
                curextremes=[calc_extreme_value(curparams,distribution,year,len(usevals)/float(nt)) for year in year_intervals]
                curparams=list(curparams)
                curparams.extend([len(usevals)/float(nt)])
                params[i/shortcut].append(curparams)
                extreme_vals[i/shortcut].append(curextremes)
            else:
                params[i/shortcut].append(list(np.array(curparams)*0))
                extreme_vals[i/shortcut].append(np.array(curextremes)*0)
            
    print("Writing Output")
    outputparams=np.array(params)
    extremes=np.array(extreme_vals)
    distname=str(distribution).split('.')[3].split()[0]
    
    if huclist!=[]:
        print(huclist.shape)
        print(extremes.shape)
        extravars=[Bunch(data=huclist,name="hucs",dims=('y',),dtype='d',attributes=None)]
    else:
        extravars=None
        
    swim_io.write(outputname+"_params"+"_"+distname,outputparams,extravars=extravars)
    swim_io.write(outputname+"_values"+"_"+distname,extremes,extravars=extravars)
    maxvals=[50,100,150,200]
    for i in range(extremes.shape[2]):
        plt.clf()
        if extremes.shape[0]==1:
            plt.imshow(extremes[:,:,i].repeat(extremes.shape[1],axis=0),
                       vmin=0,vmax=maxvals[i]*0.66,origin="lower")
        else:
            plt.imshow(extremes[:,:,i],vmin=0,vmax=maxvals[i]*0.66,origin="lower")
        plt.title(str(year_intervals[i])+"yr return interval")
        plt.colorbar()
        plt.savefig(outputname+"_yr{0:03d}.png".format(year_intervals[i]))
        
    
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Calculate nday/nyear return period events. ')
        parser.add_argument('-method',dest="methods",nargs="?",action='store',
                    default=["SDmon","CAe0","SARe0","SDe0","CA","CAe1","SDe1","SARe1","SDe","CAe","SD"])
        parser.add_argument('-model',dest="models",nargs="?",action='store',
                    default=["ncep","narr"])
        parser.add_argument('-res',dest="resolution",nargs="?",action='store',
                    default=["12km","6km"])
        parser.add_argument('-var',dest="variable",nargs="?",action='store',
                    default=["pr"])
        parser.add_argument('-bc',dest="BC",nargs="?",action='store',
                    default=["BC",""])
        parser.add_argument('-agg',dest="aggfactor",nargs="?",action='store',
                    default=0)
        parser.add_argument('-huc',dest="huc",nargs="?",action='store',
                    default="")
        parser.add_argument('-out',dest="outputdir",nargs="?",action='store',
                    default="./")
        parser.add_argument('-yearsearch',dest="yearsearch",nargs="?",action='store',
                    default=["200[1-2]"])
        parser.add_argument ('--runobs', action='store_true',
                default=False, help='verbose output', dest='runobs')
        
        parser.add_argument ('-distribution', action='store_true',
                default="gamma", help='verbose output', dest='distribution')
        
        parser.add_argument('-v', '--version',action='version',
                version='extreme_pr v.1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        if args.distribution=="gamma":
            distribution=stats.gamma
        elif args.distribution=="weibull":
            distribution=stats.weibull_min
        elif args.distribution=="exponential":
            distribution=stats.expon
        else:
            distribution=stats.gamma
        
        geosubset=[35,43,-113,-101]
        runobs=args.runobs
        runstat=not runobs
        aggfactor=args.aggfactor
        outputdir=args.outputdir
        hucfile=args.huc
        pr_threshold=0.5
        print(outputdir)
        # stats_driver.drive requires lists to iterate over, but CLI args will be individual
        #  elements.  However, default values are all lists, so we don't want to make them
        # lists of lists so we have to test to see if they are a list already first...
        for k in args.__dict__.keys():
            if type(args.__dict__[k])!=list:
                args.__dict__[k]=[args.__dict__[k]]
        
        if runobs:
            print("Running Observations Only")
            sys.stdout.flush()
        stats_driver.drive(calc_extremes,
                    yearsearch=args.yearsearch[0],
                    obs=runobs,stat=runstat,
                    extra=[distribution,pr_threshold,geosubset,aggfactor,outputdir,hucfile],
                    methods=args.methods,
                    BCs=args.BC,
                    models=args.models,
                    resolutions=args.resolution,
                    variables=args.variable,
                    extendedmethods=False)
                    
    except Exception as e:
        print("Error Unexpected Exception")
        print(e)
        traceback.print_exc()
        os._exit(1)
