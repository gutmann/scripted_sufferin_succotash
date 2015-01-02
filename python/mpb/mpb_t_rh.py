#!/usr/bin/env python

"""
SYNOPSIS

    mpb_t-rh.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION

    Extracts Temperature and Relative Humidity data from a data file
    Fills missing date/data values if possible. 
    Outputs a new file that only contains T, RH pairs of data
    
EXAMPLES

    mpt_t-rh.py grey_cat.txt grey_t-rh.txt

EXIT STATUS

    None

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    1.0

    
"""

import sys
import os
import traceback
import argparse

import numpy as np
from iterdim import iterdim

import date_fun
import load_data

def filldates(d):
    mjds=date_fun.datearr2mjd(d[:,0:5])
    dt=np.median(np.diff(mjds))
    for i in range(1,len(mjds)):
        if mjds[i]<40000:
            mjds[i]=mjds[i-1]+dt
    if mjds[0]<40000:
        mjds[0]=mjds[1]-dt
    dates=date_fun.mjd2date(mjds,roundseconds=True)
    d[:,0:5]=dates[:,0:5]

def writedata(d,filename):
    date_cols=np.array([0,1,2,3,4])
    tcols=np.array([6,8,19,21,29,31])
    rhcols=tcols-1
    datacols=np.array([5,6,7,8,18,19,20,21,28,29,30,31])
    fmt="{0[0]:5.0f} {0[1]:3.0f} {0[2]:3.0f} {0[3]:3.0f} {0[4]:3.0f}"
    header="Year Month Day Hour Minute"
    for i in range(0,len(datacols),2):
        fmt+=" {0["+str(datacols[i])+"]:8.3f} {0["+str(datacols[i+1])+"]:8.3f}"
        header+=" Rel.Humidity Temperature"
    fmt+=" {0["+str(d.shape[1])+"]:8.3f} {0["+str(d.shape[1]+1)+"]:8.3f}"
    header+=" Avg.Rel.Humidity Avg.Temperature"
    fmt+='\n'
    header+='\n'
    maskeddata=np.ma.masked_values(d,-99)
    outputd=np.zeros((d.shape[0],d.shape[1]+2))
    outputd[:,:-2]=d
    outputd[:,-2]=np.mean(maskeddata[:,rhcols],axis=1)
    outputd[:,-1]=np.mean(maskeddata[:,tcols],axis=1)
    with open(filename,'w') as f:
        f.write(header)
        for currentline in iterdim(outputd):
            f.write(fmt.format(currentline))

# def tIsGood(d):
#     return (np.abs(d).min(axis=1)>0) & (np.abs(d).max(axis=1)<40)
# 
# def rhIsGood(d):
#     return (d.min(axis=1)>=0) & (d.max(axis=1)<=1)
    
    
def filldata(d):
    tcols=np.array([6,8,19,21,29,31])
    rhcols=tcols-1
    oldd=d.copy()
    for this in iterdim(d):
        tmp=np.where((this[rhcols]<0.05) | (this[rhcols]>0.99))
        if len(tmp[0])>0:
            this[rhcols[tmp]]=-99
        tmp=np.where((this[tcols]<-40) | (this[tcols]>40))
        if len(tmp[0])>0:
            this[tcols[tmp]]=-99
    
    for var in iterdim(d[:,rhcols],axis=1):
        delta=np.diff(var)
        curve=np.abs(delta[1:]-delta[:-1])
        tmp=np.where(curve>0.1)
        if len(tmp[0])>0:
            var[tmp[0]+1]=-99
    for var in iterdim(d[:,tcols],axis=1):
        delta=np.abs(np.diff(var))
        curve=delta[1:]+delta[:-1]
        tmp=np.where(curve>10)
        if len(tmp[0])>0:
            var[tmp[0]+1]=-99
        
    print("ASSUMING specific red 2011 data!!")
    print("NEED and if --red flag (or something) then set these by date not position")
    d[25000:26850,rhcols[[3,5]]]=-99
    d[31410:31703,rhcols[4]]=-99
    d[25000:26850,tcols[[3,5]]]=-99
    d[31410:31703,tcols[4]]=-99
    mjds=date_fun.datearr2mjd(d[:,0:5])
    # d[40600:50000,rhcols[4]]=-99
    
    for i in range(1,d.shape[0]-1):
        thisline=oldd[i-1:i+2,:].sum(axis=0)
        tmp=np.where(thisline==0)
        if len(tmp[0])>0:
            d[i-1:i+2,tmp[0]]=-99
            
            
        goodt=list(tcols[np.where(d[i,tcols]>-50)[0]])
        initmean=np.mean(d[i,goodt])
        deltas=0-np.abs(d[i,goodt]-initmean)
        tmp=np.array(sorted(zip(deltas,goodt)))
        try:
            sortedgood=list(tmp[:,1])
            for r in sortedgood:
                sortedgood.remove(r)
                if len(sortedgood)>1:
                    initadj=7-len(sortedgood)
                    adjustment=np.choose(initadj>=3,(3,initadj))
                    serr=np.std(d[i,sortedgood])*adjustment
                    maxerr=np.choose(serr<2,(serr,2))
                    maxerr=np.choose(maxerr>5,(maxerr,5))
                    if np.abs(d[i,r]-np.mean(d[i,sortedgood]))>maxerr:
                        d[i,r]=-99
                    else:
                        break
                        sortedgood.append(r)
                else:
                    break
                    sortedgood.append(r)
        except:
            pass

        goodrh=list(rhcols[np.where(d[i,rhcols]>-50)[0]])
        initmean=np.mean(d[i,goodrh])
        deltas=0-np.abs(d[i,goodrh]-initmean)
        tmp=np.array(sorted(zip(deltas,goodrh)))
        try:
            sortedgood=list(tmp[:,1])
            for r in sortedgood:
                sortedgood.remove(r)
                if len(sortedgood)>1:
                    initadj=7-len(sortedgood)
                    adjustment=np.choose(initadj>=3,(3,initadj))
                    serr=np.std(d[i,sortedgood])*adjustment
                    maxerr=np.choose(serr<0.05,(serr,0.05))
                    maxerr=np.choose(maxerr>0.1,(maxerr,0.1))
                    # print(maxerr)
                    if np.abs(d[i,r]-np.mean(d[i,sortedgood]))>maxerr:
                        d[i,r]=-99
                    else:
                        break
                        sortedgood.append(r)
                else:
                    break
                    sortedgood.append(r)
        except:
            pass
            # print(tmp)
            
    # I could be doing some complicated backfilling, but for now just mark bad
    # tbad=np.where(~tIsGood(d[:,tcols]))
    # rhbad=np.where(~rhIsGood(d[:,rhcols]))
    # for t in tcols:
    #     bad=np.where(not tIsGood(d[:,t]))
    #     if len(bad[0])>0:
    #         fill(d[bad[0],tcols],t,d[:,tcols])
  

def main (filename,outputfile):
    
    d=load_data.cols(filename)
    filldates(d)
    filldata(d)
    writedata(d,outputfile)
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Pull T, RH data from a set of files (for CP only). ')
        parser.add_argument('filename',help='name of file to read')
        parser.add_argument('outputfilename',help='name of file to write')
        parser.add_argument('-v', '--version',action='version',
                version='mpb_t-rh 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        exit_code = main(args.filename,args.outputfilename)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)
