#!/usr/bin/env python
# encoding: utf-8
"""
combine_decagon.py

Created by Ethan Gutmann on 2011-09-13.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""
# std libraries
import sys
import getopt
# std third party libs
import numpy as np
import glob
# my libs
import load_data
from bunch import Bunch

help_message = '''
Combine a sequence of Decagon datalogger files

combine_decagon filesearch outputfile
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def load_all_data(files):
    alldata=list()
    n=0
    i=0
    ncols=0
    mindate=999999
    maxdate=0
    for thisfile in files:
        print(thisfile)
        alldata.append(load_data.decagon(thisfile))
        n+=alldata[i].shape[0]
        thiscols=alldata[i].shape[1]
        curmin=np.min(alldata[i][:,0])
        curmax=np.max(alldata[i][:,0])
        mindate=np.choose(mindate>curmin,(mindate,curmin))
        maxdate=np.choose(maxdate<curmax,(maxdate,curmax))
        ncols=np.choose(ncols<thiscols,(ncols,thiscols))
        i+=1
    return Bunch(mindate=mindate,maxdate=maxdate,ncols=ncols,data=alldata,n=n)


def write_alldata(data,filename):
    sz=data.shape
    with open(filename,'wt') as f:
        for i in range(sz[0]):
            for j in range(sz[1]):
                f.write(str(data[i,j])+' ')
            f.write('\n')


def combine_decagon(filesearch,outputfile):
    files=glob.glob(filesearch+'*')
    data=load_all_data(files)

    dt=(data.data[0][-1,0]-data.data[0][0,0])/(data.data[0].shape[0]-1)
    n=np.round((data.maxdate-data.mindate)/dt)
    outputdata=np.zeros((n+1,data.ncols))
    for i in range(len(files)):
        curlen=data.data[i].shape[0]
        ncols=data.data[i].shape[1]
        # unforuntately, data aren't always evenly spaced, 
        # so we have to compute each point independantly
        # if points are evenly spaced you could comment out the for loop
        # and use the commented section below the loop
        for j in range(curlen):
            curpt=np.round((data.data[i][j,0]-data.mindate)/dt)
            outputdata[curpt,:ncols]=data.data[i][j,:]
        # startpt=np.round((data.data[i][0,0]-data.mindate)/dt)
        # endpt=np.round((data.data[i][-1,0]-data.mindate)/dt)
        # outputdata[startpt:endpt+1,:ncols]=data.data[i]
        
    write_alldata(outputdata,outputfile)
    


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
        if len(args)<2:
            raise Usage(help_message)
            
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    combine_decagon(args[0],args[1])

if __name__ == "__main__":
    sys.exit(main())
