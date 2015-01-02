#!/usr/bin/env python
# encoding: utf-8
"""
ingest_data.py

Created by Ethan Gutmann on 2011-09-14.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""

import sys
import getopt
import subprocess
import re
import os
import glob

import numpy as np

import load_data
import date_fun
from bunch import Bunch

help_message = '''
combine one or more (initially only chimney park campbell sci) datalogger file with a master file

ingest_data.py <filesearchpattern> [masterfilename <default="master">]

filesearchpattern can be a single filename
output master file will actually be 3 or 4 files with _data10x.txt suffixes
'''
defaultmasterfile='master'
filesuffix=['_data101.txt','_data102.txt','_data104.txt','_data105.txt']

def load_masterdata(masterfile):
    # masterdata files should all be simple load_data.cols-able
    masterdata101=load_data.cols(masterfile+filesuffix[0])
    masterdata102=load_data.cols(masterfile+filesuffix[1])
    masterdata104=load_data.cols(masterfile+filesuffix[2])
    try: 
        masterdata105=load_data.cols(masterfile+filesuffix[3])
    except IOError:
        print("Missing data from tbl105 in file: "+masterfile)
        masterdata105=None
    
    return Bunch(data101=masterdata101,data102=masterdata102,
                data104=masterdata104,data105=masterdata105)

# this is one important piece of work combing the two datasets
def combine_data(newdata,olddata):
    # read modified julian dates for each
    newmjd=date_fun.date2mjd(newdata[:,0],newdata[:,1],newdata[:,2],newdata[:,3],newdata[:,4])
    oldmjd=date_fun.date2mjd(olddata[:,0],olddata[:,1],olddata[:,2],olddata[:,3],olddata[:,4])
    
    tmp=np.where(newmjd>50000)
    if tmp[0].size>1:
        newmjd=newmjd[tmp]
        newdata=newdata[tmp[0],:]
    
    # find the bounds on the output data set in date space
    mindate=min((min(newmjd),oldmjd[0]))
    mindate=max((mindate,date_fun.date2mjd(2007,1,1,0,0)))
    maxdate=max((max(newmjd),oldmjd[-1]))
    maxdate=min((maxdate,date_fun.date2mjd(2013,1,1,0,0)))
    # find dt using old mjd because it is "guaranteed" to be uniformly spaced (no missing lines)
    # ... though now I think newmjd should be too (by make_even_dates())
    # dt=(oldmjd[-1]-oldmjd[0])/(oldmjd.size-1)
    dt=np.double(np.median(oldmjd[1:]-oldmjd[:-1]))
    dt=np.double(15.0)/24.0/60.0
    # convert bounds to length of data
    datalen=np.round((maxdate-mindate)/dt)
    # create a new array for the output dataset filling bounds
    finaldata=np.zeros((datalen+1,olddata.shape[1]))
    # find the old data start and end points so we can plunk the old data into the output data
    oldstart=np.round((oldmjd[0]-mindate)/dt)
    oldend=np.round((oldmjd[-1]-mindate)/dt)
    finaldata[oldstart:oldend+1,:]=olddata
    # then loop through the new data putting it in the proper place, 
    # assumes nothing about the order or continuity of the new data 
    # can be done in two steps: 
    points=np.round((newmjd-mindate)/dt)
    finaldata[points.astype('i'),:]=newdata
    # fill in new, potentially missing, dates
    outputmjd=np.arange(datalen+1)*dt+mindate
    outputdates=date_fun.mjd2date(np.double(np.round(outputmjd/dt))*dt,roundseconds=True)
    finaldata[:,0:5]=outputdates[:,0:5]
    # for i in range(newdata.shape[0]):
    #   curpt=(newmjd[i]-mindate)/dt
    #   if (curpt>0) & (curpt<finaldata.shape[0]):
    #       finaldata[curpt,:]=newdata[i,:]
    
    return finaldata


# no longer used after discovering that some files arbitrarily miss a tbl or two 
# so you have to search through every line!
def greptofiles(filename):
    p1=subprocess.Popen("grep '^101' "+filename+" >"+filename+'tmp'+filesuffix[0],shell=True)
    p2=subprocess.Popen("grep '^102' "+filename+" >"+filename+'tmp'+filesuffix[1],shell=True)
    p3=subprocess.Popen("grep '^104' "+filename+" >"+filename+'tmp'+filesuffix[2],shell=True)
    p4=subprocess.Popen("grep '^105' "+filename+" >"+filename+'tmp'+filesuffix[3],shell=True)
    p1.wait()
    p1=subprocess.Popen("comma2space "+filename+'tmp'+filesuffix[0]+' '+filename+filesuffix[0],shell=True)
    p2.wait()
    p2=subprocess.Popen("comma2space "+filename+'tmp'+filesuffix[1]+' '+filename+filesuffix[1],shell=True)
    p3.wait()
    p3=subprocess.Popen("comma2space "+filename+'tmp'+filesuffix[2]+' '+filename+filesuffix[2],shell=True)
    p4.wait()
    p4=subprocess.Popen("comma2space "+filename+'tmp'+filesuffix[3]+' '+filename+filesuffix[3],shell=True)
    p1.wait()
    os.remove(filename+'tmp'+filesuffix[0])
    p2.wait()
    os.remove(filename+'tmp'+filesuffix[1])
    p3.wait()
    os.remove(filename+'tmp'+filesuffix[2])
    p4.wait()
    os.remove(filename+'tmp'+filesuffix[3])


def load_raw_data(filename):
    cursize=10000
    mindate=date_fun.ydoyhm2mjd(2008, 300, 2400)
    # could probably now be converted to a with open(filename,'rt') as f: but won't work with python<2.5... so?
    f=open(filename,'rt')
    try:
        line=f.next()
        while (re.match('^101,',line)==None):
            line=f.next()
        curdata=np.array(line.split(',')).astype('f')
        mjd=date_fun.ydoyhm2mjd(curdata[1],curdata[2],curdata[3])
        
        if mjd<mindate:
            mjd+=1102.71875# 1102.7194444
        dt=np.double(15.0)/60/24.0# np.median(mjd[1:]-mjd[:-1])
        mjd=np.round(mjd/dt)*dt
            
        date=(date_fun.mjd2date(mjd,roundseconds=True))[:-1]
        data101=np.zeros((cursize,curdata.size+1))
        data101[0,:]=np.hstack([date,curdata[4:]])
        inputdata=[data101,None,None,None,None]
        state=0
        i=0
        statestr=['^101','^102','^103','^104','^105']
        # search through the file, at EOF we will break with a StopIteration exception
        while True:
            line=f.next()
            # update the current state (i.e. which table 101,102,... we should be in)
            state=(state+1) % 5
            # we we are indeed in the right table, process this line
            if (re.match(statestr[state],line)!=None):
                # if this line is a tbl 101 (includes a date) process it separately
                if state==0:
                    # we advance to the next data group (i+=1)
                    i+=1
                    # if we spilled over the current max size, double the size of all groups
                    if i>=cursize:
                        for j in range(len(statestr)):
                            tmp=inputdata[j].copy()
                            inputdata[j]=np.zeros((cursize*2,tmp.shape[1]))
                            inputdata[j][0:cursize,:]=tmp
                        cursize*=2
                    # read the data as a line of floats
                    curdata=np.array(line.split(',')).astype('f')
                    # read the date into a modified julian day
                    mjd=date_fun.ydoyhm2mjd(curdata[1],curdata[2],curdata[3])
                    # correct for a possible date error in some of the 2011 data
                    if mjd<mindate:
                        mjd+=1102.71875# 1102.7194444
                    # store the date as Year,Month,Day,Hour,Minute
                    # convert to even dt intervals
                    # dt=np.double(15.0)/60/24.0# np.median(mjd[1:]-mjd[:-1])
                    mjd=np.round(mjd/dt)*dt
                    
                    date=(date_fun.mjd2date(mjd,roundseconds=True))[:-1]
                    # stack the new date and the existing data into the current input data point
                    inputdata[state][i,:]=np.hstack([date,curdata[4:]])
                
                else:
                    # if we are in any of the other tables process them "normally"
                    curdata=np.array(line.split(',')).astype('f')
                    # if we have never entered this tbl before create an array for it
                    if inputdata[state]==None:
                        inputdata[state]=np.zeros((cursize,curdata.size+4))
                    # stack the date and the current data into the right spot in the input data
                    inputdata[state][i,:]=np.hstack((date,curdata[1:]))
            # if we did not match the correct table for this line, search through lines until we hit
            # a tbl 101 line again so we know we have the correct date. 
            else:
                # first print the offending line so I can look at the input file
                if (inputdata[state]!=None): print(line)
                # now loop until we match tbl 101
                while (re.match('^101,',line)==None):
                    line=f.next()
                # we found a date, so advance the date group pointer i
                i+=1
                # if we spilled over the current max size, double the size of all groups
                if i>=cursize:
                    for j in range(len(statestr)):
                        tmp=inputdata[j].copy()
                        inputdata[j]=np.zeros((cursize*2,tmp.shape[1]))
                        inputdata[j][0:cursize,:]=tmp
                    cursize*=2
                # and our current state is tbl 101
                state=0
                # then process normally (store data, get and correct date, store date+data)
                curdata=np.array(line.split(',')).astype('f')
                # get date
                mjd=date_fun.ydoyhm2mjd(curdata[1],curdata[2],curdata[3])
                # correct date if necessary
                if mjd<mindate:
                    mjd+=1102.71875
                mjd=np.round(mjd/dt)*dt
                # convert modified julian day to Year,Month,Day,Hour,Minute
                date=(date_fun.mjd2date(mjd,roundseconds=True))[:-1]
                # stack date and data into input data
                inputdata[state][i,:]=np.hstack([date,curdata[4:]])
                
                
                
    except StopIteration:
        f.close()
    
    if inputdata:
        return Bunch(data101=inputdata[0],data102=inputdata[1],
                data104=inputdata[3],data105=inputdata[4])
    else:
        return None
    
        


# take a potentially irregularly spaced (in time) dataset and add missing points if necessary
def make_even_dates(data,mjdinput=None):
    # make dates into a single number we can work with (modified julian day)
    if mjdinput==None:
        mjd=date_fun.date2mjd(data[:,0],data[:,1],data[:,2],data[:,3],data[:,4])
    else:
        mjd=mjdinput
    # subset down to periods where the date is valid (post ~1995 in this case)
    tmp=np.where(mjd>50000)
    if tmp[0].size>1:
        mjd=mjd[tmp]
        data=data[tmp[0],:]
    # find the min and max times in the dataset
    mindate=min(mjd)
    mindate=max((mindate,date_fun.date2mjd(2007,1,1,0,0)))
    maxdate=max(mjd)
    maxdate=min((maxdate,date_fun.date2mjd(2013,1,1,0,0)))
    # find the typical (median) time step
    dt=np.double(np.round(np.median(mjd[1:]-mjd[:-1])*24*60))/24.0/60
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
    if mjdinput==None:
        outputdates=date_fun.mjd2date(outputmjd,roundseconds=True)
        outputdata[:,0:5]=outputdates[:,0:5]
    else:
        outputdata[:,0]=outputmjd
    return outputdata
    
    

# this is one large piece of work take an awkward file format and parse it into 4 reasonable formats
def load_newdata(filename):
    
    data=load_raw_data(filename)
    if data==None:
        return None
    for key in data.keys():
        if data[key]!=None:
            data[key]=make_even_dates(data[key])
    return data
    
    # greptofiles(filename) 
    # data101=data.data101
    # # data101=load_data.cols(filename+filesuffix[0])
    # mjd=date_fun.ydoyhm2mjd(data101[:,1],data101[:,2],data101[:,3])
    # 
    # # the red / upper / _2 site had the clock off by 1102.72 days for a period in 2011
    # # this just fixes that offset
    # tmp=np.where(mjd<date_fun.ydoyhm2mjd(2008, 300, 2400))
    # if tmp[0].size>0:
    #   mjd[tmp]+=1102.7194444
    # 
    # #roundseconds prevents it from coming up with seconds=59 and minutes=minutes-1
    # dates=date_fun.mjd2date(mjd,roundseconds=True) 
    # dates=dates[:,:-1] #lopp off the seconds
    # 
    # data101=np.hstack((dates,data.data101[:,4:]))
    # data102=np.hstack((dates,data.data102[:,1:]))
    # data104=np.hstack((dates,data.data104[:,1:]))
    # if data.data105!=None:
    #   data105=np.hstack((dates,data.data105[:,1:]))
    # else:
    #   data105=None
    # 
    # return Bunch(data101=data101,data102=data102,
    #           data104=data104,data105=data105)

# loop through data.keys writing a file for each key
# loop through i,j in each dataset writing data to file
def write_master_data(filename,data):
    for key in data.keys():
        curdata=data[key]
        if curdata!=None:
            fileoutput=open(filename+'_'+key+'.txt','wt')
            sz=curdata.shape
            for i in range(sz[0]):
                for j in range(sz[1]):
                    fileoutput.write(str(curdata[i,j])+' ')
                fileoutput.write('\n')
            fileoutput.close()


def mpb_ingest_data(filename,masterfile=None):
    files=glob.glob(filename)
    # we can import more than one file at once if we were given a pattern that matches multiple files
    if len(files)>1:
        for f in files:
            print(f)
            mpb_ingest_data(f,masterfile=masterfile)
            if masterfile==None:masterfile=defaultmasterfile
        return
    # in case we were given a pattern that only matches one file
    else: 
        filename=files[0]
    # this is what ingest data really does
    # load data into a dict with keys for each table
    newdata=load_newdata(filename)
    if newdata==None:
        return
    # if we were given a master filename try to process it
    if masterfile!=None:
        # if we were given a master file, but it/they don't exist, just make masterdata=newdata
        masterfiles=glob.glob(masterfile+'*')
        if len(masterfiles)>0:
            # if we do have masterfiles, then read masterfiles and combine them with newdata
            # loads master data into a dict with keys for each table
            masterdata=load_masterdata(masterfile)
            for key in masterdata.keys():
                if masterdata[key]==None:
                    masterdata[key]=newdata[key]
                elif newdata[key]==None:
                    pass
                else:
                    # this is the main ingestion/combination routine
                    masterdata[key]=combine_data(newdata[key],masterdata[key])
        else:
            masterdata=newdata
    else:
        masterdata=newdata
    if masterfile==None:
        masterfile=defaultmasterfile
    # finally write the new master file
    write_master_data(masterfile,masterdata)


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
    
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
        if len(args)==0:
            raise Usage(help_message)
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    if len(args)==2:
        mpb_ingest_data(args[0],masterfile=args[1])
    if len(args)==1:
        mpb_ingest_data(args[0])


if __name__ == "__main__":
    sys.exit(main())
