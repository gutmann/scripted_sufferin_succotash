#!/usr/bin/env python

"""
SYNOPSIS

    wrf2daily.py [-h] [--verbose] [-v, --version] <filename_search>

DESCRIPTION

    Searches for a set of WRF output files and converts precip and temperture
    to daily values (daily total precip and daily min and max temperatures.)
    
EXAMPLES

    wrf2daily.py 'wrf\*.nc'

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
# import argparse
import glob

import numpy as np

import swim_io
from bunch import Bunch

UTC_offset=7
date_var="Times"
data_var="T2"

def make_output(filename,nfiles):
    d=swim_io.read_nc(filename,var=data_var,returnNCvar=True)
    sz=d.data.shape
    d.ncfile.close()
    nsteps=int(np.round(sz[0]/24.0))
    outputdata=np.zeros((nfiles*nsteps,2,sz[-2],sz[-1]))
    print(nsteps,outputdata.shape)
    return Bunch(dates=[],
                data=outputdata,
                nstepsperfile=nsteps)

def get_data_for_var(file1,file2,offset,var="T2"):
    d=swim_io.read_nc(file1,var=var,returnNCvar=True)
    if d.data.shape[0]>offset+24:
        output=d.data[offset:offset+24,...]
        d.ncfile.close()
        return (date,output)
    else:
        d2=swim_io.read_nc(file2,var=var,returnNCvar=True)
        output=np.vstack([d.data[offset:,...],d2.data[:24-offset,...]])
        d.ncfile.close()
        d2.ncfile.close()
        return output

def dailyfromfiles(file1,file2,output,time_offset=0,curpos=0):
    for i in range(output.nstepsperfile):
        temperature=get_data_for_var(file1,file2,time_offset+(i*24),var=data_var)
        output.data[curpos,0,:,:]=np.min(temperature,axis=0)
        output.data[curpos,1,:,:]=np.max(temperature,axis=0)
        curdate=''.join(swim_io.read_nc(file1,var=date_var).data[i*24])[0:10].replace('-','')
        output.dates.append(long(curdate))
        curpos+=1
    return curpos
    
def write_output(outputfilename,outputdata):
    sz=outputdata.data.shape
    ncout=swim_io.NC_writer(outputfilename+"_tasmin",sz[3],sz[2],var="tasmin")
    for i in range(sz[0]):
        ncout.appendToVar(outputdata.data[i,0,...],date=outputdata.dates[i])
    ncout.close()
    ncout=swim_io.NC_writer(outputfilename+"_tasmax",sz[3],sz[2],var="tasmax")
    for i in range(sz[0]):
        ncout.appendToVar(outputdata.data[i,1,...],date=outputdata.dates[i])
    ncout.close()

def main (filename):
    UTC_offset=7
    files=np.sort(glob.glob(filename))
    outputdata=make_output(files[0],len(files))
    cur_position=0
    for f1,f2 in zip(files[:-1],files[1:]):
        print(f1)
        cur_position=dailyfromfiles(f1,f2,outputdata,
                time_offset=UTC_offset,curpos=cur_position)
    outputdata.data=outputdata.data[:cur_position,...]
    write_output("/glade/scratch/gutmann/wrf_daily",outputdata)
    
    
    
if __name__ == '__main__':
    main("wrfout*.nc")
    # try:
    #     parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ')
    #     parser.add_argument('filename',action='store')
    #     parser.add_argument('-v', '--version',action='version',
    #             version='0.1')
    #     parser.add_argument ('--verbose', action='store_true',
    #             default=False, help='verbose output', dest='verbose')
    #     args = parser.parse_args()
    # 
    #     exit_code = main(args.filename)
    #     if exit_code is None:
    #         exit_code = 0
    #     sys.exit(exit_code)
    # except KeyboardInterrupt, e: # Ctrl-C
    #     raise e
    # except SystemExit, e: # sys.exit()
    #     raise e
    # except Exception, e:
    #     print('ERROR, UNEXPECTED EXCEPTION')
    #     print(str(e))
    #     traceback.print_exc()
    #     os._exit(1)
