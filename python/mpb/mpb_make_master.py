#!/usr/bin/env python
import sys, subprocess, os

import numpy as np
from iterdim import iterdim

import load_data, date_fun

def convertcommas2spaces(filename):
    outputfile=filename+'.txt'
    p=subprocess.Popen('comma2space '+filename+' '+outputfile,shell=True)
    p.wait()
    return outputfile

def fetch_mjd(data):
    year=data[:,1]
    doy=data[:,2]
    hhmm=data[:,3]
    dates=date_fun.ydoyhm2mjd(year,doy,hhmm)
    return dates

# take a potentially irregularly spaced (in time) dataset and add missing points if necessary
def make_even_dates(mjd,data):
    
    # subset down to periods where the date is valid (post ~1995 in this case)
    tmp=np.where(mjd>50000)
    if tmp[0].size>1:
        mjd=mjd[tmp]
        data=data[tmp[0],:]
    # find the min and max times in the dataset
    mindate=np.min(mjd)
    maxdate=np.max(mjd)
    # find the typical (median) time step
    # dt=np.double(np.round(np.median(mjd[1:]-mjd[:-1])*24*60))/24.0/60
    dt=np.double(15)/60/24.0
    # create a new output data set that spans from mindate to max date in dt time steps
    outputdata=np.zeros((np.round((maxdate-mindate)/dt+1),data.shape[1]))
    # compute where each point in the raw dataset falls in the output dataset
    points=np.round((mjd-mindate)/dt)
    # put the raw date in its proper place in the output dataset
    outputdata[points.astype('i'),:]=data
    # now fill in all of the dates (so we don't have 0 0 0 0 0 for some dates)
    # first compute a linearly increases range of dates
    outputmjd=np.arange(mindate,maxdate+dt/10.0,dt).astype(np.double)
    outputmjd=np.round(outputmjd/dt)*dt
    # convert those MJDs to an array of dates and put them into the output dataset
    outputdates=date_fun.mjd2date(outputmjd,roundseconds=True)
    outputdata[:,0:5]=outputdates[:,0:5]
    return outputdata

def write_masterfile(filename,data):
    sz=data.shape
    with open(filename,'wt') as f:
        for i in range(sz[0]):
            for j in range(sz[1]):
                f.write(str(data[i,j])+' ')
            f.write('\n')

def make_master(filename, masterfile=None):
    
    if masterfile==None:
        masterfile='master'+filename+'.txt'
    spacefile=convertcommas2spaces(filename)
    data=load_data.cols(spacefile)
    os.remove(spacefile)
    
    mjds=fetch_mjd(data)
    data=make_even_dates(mjds,data)
    # np.savetxt(master,data)
    write_masterfile(masterfile,data[:,:9])

if len(sys.argv)==2:
    make_master(sys.argv[1])
elif len(sys.argv)==3:
    make_master(sys.argv[1],sys.argv[2])
else:
    print("make_master.py <filename> [masterfilename]")
