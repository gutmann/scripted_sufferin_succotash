#!/usr/bin/env python
"""
SYNOPSIS

    stat_histograms.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION

    Takes an input file that describes the points/files/variables to read
    Generates histograms for each point or set of points for each variable/file
    Outputs data and plots of said histograms. 
    
EXAMPLES

    stat_histograms.py histo.info

EXIT STATUS

    None

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    0.1
    
"""

import sys
import os
import traceback
import argparse
import glob
import re

import numpy as np
import matplotlib.pyplot as plt

import swim_io as mygis
import swim_io as nc
from bunch import Bunch

def make_info(inputline):
    '''Parse one line of input into a structure'''
    searchpath="*.nc"
    varname="pr"
    latvar="lat"
    lonvar="lon"
    minval=0.1
    maxval=100.0
    islog=1
    geofile=""
    varname="pr"
    outputname="test"
    info=Bunch(search=searchpath,latvar=latvar,lonvar=lonvar,outputname=outputname,
                minval=minval,maxval=maxval,islog=islog,geofile=geofile,varname=varname)
    for i in inputline:
        key,value=i.split("=")
        key=key.strip()
        value=value.strip()
        info[key]=type(info[key])(value)
    info.islog=(info.islog==1)
    return info

def get_continuous_line(f):
    '''Read lines from f as long as there is an & at the end of the line'''
    l=f.next().strip()
    newl=l
    while newl[-1]=="&":
        l=l[:-1]
        newl=f.next().strip()
        l+=newl
    return l

def parse_inputfile(filename):
    '''Parse inputfile into lat,lon and a series of file info structures
    
    Note all data are comma separated.  
    Within a comma separated expression (e.g. between two commas)
        ":" indicates that a range between two values should be used with a step size 
            described after a second ":"
        ";" are used for a list of values
        
    lines ending with an & (white space is ignored) are continued on the next line
    "info" structure is a series of comma separated key=value pairs. '''
    f=open(filename,'r')
    # split the lat and lon lines into comma separated groups
    lat=get_continuous_line(f).split(",")
    lon=get_continuous_line(f).split(",")
    outputlat=list()
    outputlon=list()
    for lat1,lon1 in zip(lat,lon):
        # within a comma delimited group look for ranges (: separated)
        if re.match(".*:.*",lat1):
            latinfo=np.array(lat1.split(":"),dtype="f")
            loninfo=np.array(lon1.split(":"),dtype="f")
            lat_1d=np.arange(*latinfo)
            lon_1d=np.arange(*loninfo)
            lat2d,lon2d=np.meshgrid(lat_1d,lon_1d)
            lat2d.shape=(lat2d.size)
            lon2d.shape=(lon2d.size)
            outputlat.append(lat2d)
            outputlon.append(lon2d)
        else:
            # otherwise, check for ; separated list
            curlat=lat1.split(";")
            curlon=lon1.split(";")
            if len(curlat)>1:
                lat1d=np.array(curlat,dtype="f")
                lon1d=np.array(curlon,dtype="f")
                outputlat.append(lat1d)
                outputlon.append(lon1d)
            else:
                # otherwise it should just be a single value
                outputlat.append(float(lat1))
                outputlon.append(float(lon1))
    
    finfo=list()
    for l in f:
        l=l.strip()
        if l[-1]=="&":
            l=l[:-1]+get_continuous_line(f)
        finfo.append(make_info(l.split(",")))
    f.close()
    return Bunch(info=finfo,lat=outputlat,lon=outputlon)

def find_xy(filename,latpt,lonpt,latvar,lonvar):
    '''Find the x,y point in file corresponding the latpt,lonpt'''
    lat=mygis.read_nc(filename,var=latvar).data
    lon=mygis.read_nc(filename,var=lonvar).data
    if lon.max()>180:
        lon[lon>180]=lon[lon>180]-360
    if len(lat.shape)==1:
        if len(latpt)==1:
            x=[np.argmin(np.abs(lon-lonpt))]
            y=[np.argmin(np.abs(lat-latpt))]
        else:
            x=[np.argmin(np.abs(lon-thislonpt)) for thislonpt in lonpt]
            y=[np.argmin(np.abs(lat-thislatpt)) for thislatpt in latpt]
            
    else:
        if len(lat.shape)==3:
            lat=lat[0,:,:]
            lon=lon[0,:,:]
        if len(latpt)==1:
            dists=np.abs((lon-lonpt)**2 + (lat-latpt)**2)
            y,x=np.unravel_index(dists.argmin(),dists.shape)
            y=[y]
            x=[x]
        else:
            y=list()
            x=list()
            for thislonpt,thislatpt in zip(lonpt,latpt):
                dists=np.abs((lon-thislonpt)**2 + (lat-thislatpt)**2)
                cury,curx=np.unravel_index(dists.argmin(),dists.shape)
                x.append(curx)
                y.append(cury)
    return x,y
    

def find_all_xy(info):
    '''loop through files and lat,lon pairs finding all corresponding x,y points in files'''
    allpts=list()
    latpt=info.lat
    lonpt=info.lon
    for thisinfo in info.info:
    # thisinfo=info.info[0]
        xypts=list()
        if thisinfo.geofile:
            thisfile=thisinfo.geofile
        else:
            thisfile=glob.glob(thisinfo.search)[0]
        print(thisfile)
        for curlat,curlon in zip(latpt,lonpt):
            xypts.append(find_xy(thisfile,curlat,curlon,thisinfo.latvar,thisinfo.lonvar))
            print(thisfile,xypts[-1][1][0],xypts[-1][0][0],curlat[0],curlon[0])
        allpts.append(xypts)
    return allpts
    
def read_times(filename):
    return np.arange(10)
    
def find_times(filename,month):
    try:
        nmonths=len(month)
    except TypeError:
        month=[month]
    time_data=read_times(filename)
    use_times=np.zeros(len(time_data))
    n=0
    for i,t in enumerate(time_data):
        if t in month:
            use_times[n]=i
            n+=1
    return use_times[:n]

def calc_histogram(filesearch,varname,x,y,minval,maxval,islog=False,months=None):
    '''calculate the histogram for point x,y of varname in matching files'''
    files=glob.glob(filesearch)
    fulldata=[]
    for f in files:
        d=mygis.read_nc(f,varname,returnNCvar=True)
        if months==None:
            cur=d.data[:,y.min():y.max()+1,x.min():x.max()+1]
            fulldata.extend(cur[:,y-y.min(),x-x.min()])
        else:
            times=find_times(f,months)
            fulldata.extend(d.data[times,y,x])
        d.ncfile.close()
    fulldata=np.array(fulldata)
    if re.match(".*wrf.*_tasm.*",files[0]):
        fulldata=fulldata[:2192,...]-273.15
    fulldata=fulldata[np.where(np.isfinite(fulldata))]
    if islog:
        fulldata=fulldata[np.where(fulldata>0)]
        nbins=20
        delta=minval*0.001
        fulldata[fulldata<minval]=minval+delta
        fulldata[fulldata>maxval]=maxval-delta
        y,x=np.histogram(np.log10(fulldata),range=[np.log10(minval),np.log10(maxval)],bins=nbins)
        x=10.**x
    else:
        nbins=30
        delta=(maxval-minval)/(nbins*10)
        fulldata[fulldata<minval]=minval+delta
        fulldata[fulldata>maxval]=maxval-delta
        y,x=np.histogram(fulldata,range=[minval,maxval],bins=nbins)
    return x,y

def read_datasets(info,xypts,writefast=False):
    '''return list of histogram datasets for all xypoints in files from info list'''
    allhists=[]
    for (curxy,thisinfo) in zip(xypts,info.info):
        for (x,y),lat,lon in zip(curxy,info.lat,info.lon):
            fullhist=None
            # print(glob.glob(thisinfo.search)[0],y[0],x[0],lat[0],lon[0])
            for x1,y1 in zip(x,y):
                xvals,curhist=calc_histogram(thisinfo.search,thisinfo.varname,np.array(x1),np.array(y1),
                                            thisinfo.minval,thisinfo.maxval,islog=thisinfo.islog)
                if fullhist==None:
                    fullhist=curhist
                else:
                    fullhist+=curhist
            print("   "+thisinfo.outputname)# +"%6.1f"%np.mean(x))
            if not writefast:
                allhists.append((thisinfo.islog,thisinfo.varname,thisinfo.outputname,
                                ("%7.1f_%7.1f"%(np.mean(lon),np.mean(lat))).replace(" ",""),xvals,fullhist))
            else:
                xyname=("%7.1f_%7.1f"%(np.mean(lon),np.mean(lat))).replace(" ","")
                outname=thisinfo.outputname+"_"+xyname+".nc"
                nc.write(outname,np.vstack((xvals,np.hstack(([0],fullhist)) )),varname="hist")
    return allhists      
        
def plot_hist(x,y,logscale=False,xlabel=None):
    '''Plot the histogram described by x,y values'''
    plt.clf();
    ax=plt.subplot(111)
    if logscale:
        ax.set_xscale("log")
        ax.set_yscale("log")
    x1=x[:,np.newaxis].repeat(2,axis=1).reshape(len(x)*2)
    y1=[0]
    y1.extend(y[:,np.newaxis].repeat(2,axis=1).reshape(len(y)*2))
    y1.append(0)
    plt.plot(x1,y1)
    plt.xlim(x1.min(),x1.max())
    ax.yaxis.set_label_text("N")
    if xlabel:
        ax.xaxis.set_label_text(xlabel)
    plt.draw()
    
    
        
def make_plots(data):
    '''Plot histograms and write data files for each data set'''
    for (islog,varname,outname,xyname,x,y) in data:
        # plot_hist(x,y,logscale=islog,xlabel=varname)
        # plt.savefig(outname+"_"+xyname+".png")
        nc.write(outname+"_"+xyname,np.vstack((x,np.hstack(([0],y)) )),varname="hist")
        

def main (filename):
    writefast=True
    print("Reading Input File")
    info=parse_inputfile(filename)
    print("Calculationg positions")
    xypts=find_all_xy(info)
    print("Reading data")
    data=read_datasets(info,xypts,writefast)
    if not writefast:
        print("Writing data")
        make_plots(data)
    
    
if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Create histograms for points in data files. ')
        parser.add_argument('filename',action='store')
        parser.add_argument('-v', '--version',action='version',
                version='stat_histograms 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        exit_code = main(args.filename)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
