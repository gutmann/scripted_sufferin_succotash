#!/usr/bin/env python

"""
SYNOPSIS

    nc_addvars.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION

    TODO This describes how to use this script.
    This docstring will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION

    
"""

import sys
import os
import traceback
import argparse
import glob

import numpy as np
import swim_io
import nio
# from netCDF4 import Dataset

def add_vars_to_file(filename,inputvars,varnames,output_prefix="added_"):
    '''Add the list of data (inputvars) and names (varnames) to file (filename)
    
    Note: the dimensions of the inputvars must match some existing 
          dimension in the input file. 
    
    filename is the name of an existing netcdf file
    input vars is a list of numpy arrays
    varnames is a list of strings 
    output_prefix = prefix to be added to the new filename
            default is "added_"

    A new file will be created titled: output_prefix+filename
    '''
    outputfilename=output_prefix+filename.split("/")[-1]
    outputfilename=outputfilename.strip(".nc") #strip off the .nc term first
    # if there wasn't a directory specified, make sure we aren't specifying "/"
    if len(filename.split("/"))>1:
        outputdir="/".join(filename.split("/")[:-1])+"/"
    else:
        outputdir=""
    # if there wasn't a directory specified, make sure we aren't specifying "/"
    if len(outputdir)==1:
        outputdir=""
    outputf=nio.open_file(outputdir+outputfilename, 'w',format="nc")
    inputf=nio.open_file(filename, 'r')
    maxlen=0
    # loop over the input variables finding the maximum number of 
    # dimensions we need to know about
    for d in inputvars:
        maxlen=max(maxlen,len(d.shape))
    # create an array to hold the references to the dimensions corresponding
    #  to the dimensions of each input variable
    outputdims=np.array((len(inputvars),maxlen),dtype="S16")
    if len(outputdims.shape)==1:
        outputdims=outputdims.reshape((outputdims.shape[0],1))
    dimnum=0
    # loop over the dimensions in the input file
    for i,dim in enumerate(inputf.dimensions.keys()):
        dimlen=int(inputf.dimensions[dim]) #in netCDF4python parlence we need a len() around dimension (?)
        
        # create that dimension in the outputfile too
        outputf.create_dimension(dim,dimlen)
        # find any input variables that match this dimension
        for j,thisvar in enumerate(inputvars):
            finddim=np.where(np.array(thisvar.shape)==dimlen)
            if len(finddim[0])>0:
                # if we found a match
                # save the dimension in the output dims array
                for k in finddim[0]:
                    outputdims[j,k]=dim

    # loop through the existing variables in the input file
    # copying them to the output file (originally allowed subsetting in time)
    for v in inputf.variables:
        thisvar=inputf.variables[v]
        # create the variable
        outputf.create_variable(str(v),thisvar.typecode(),thisvar.dimensions)
        # copy the data over
        outputf.variables[str(v)][:]=inputf.variables[v][:]
            
    i=0
    # finally copy the new variables into the file as well
    for v,name in zip(inputvars,varnames):
        # create the variable
        outputf.create_variable(name,"f",tuple(outputdims[i,:len(v.shape)]))
        # store the data
        outputf.variables[name][:]=v.astype("f")
        i+=1
    
    inputf.close()
    outputf.close()
    

def main (filesearch="*.nc",copy_from_file="master.nc",vars2copy=["lat","lon"],
          ranges=None,inputvars=None,verbose=False,output_prefix="added_"):
    if ranges==None:
        ranges=[]
        for f in vars2copy:
            ranges.append([])
    
    files=glob.glob(filesearch)
    files=np.sort(files)
    
    # this isn't the most general way to do this, but I'm not likely to need an arbitrary number of dimensions
    if inputvars==None:
        inputvars=[]
        for v,minmax in zip(vars2copy,ranges):
            # if it is only a one dimensional array
            if len(minmax)==2:
                xmin=minmax[0]
                xmax=minmax[1]
                print(v,xmin,xmax)
                inputvars.append(swim_io.read_nc(copy_from_file,var=v).data[xmin:xmax])
            # if it is a two dimensional array
            elif len(minmax)==4:
                xmin=minmax[0]
                xmax=minmax[1]
                ymin=minmax[2]
                ymax=minmax[3]
                inputvars.append(swim_io.read_nc(copy_from_file,var=v).data[ymin:ymax,xmin:xmax])
            # three dimensions
            elif len(minmax)==6:
                xmin=minmax[0]
                xmax=minmax[1]
                ymin=minmax[2]
                ymax=minmax[3]
                tmin=minmax[4]
                tmax=minmax[5]
                inputvars.append(swim_io.read_nc(copy_from_file,var=v).data[tmin:tmax,ymin:ymax,xmin:xmax])
            # four dimensions (e.g. x,y,z,time)
            elif len(minmax)==8:
                xmin=minmax[0]
                xmax=minmax[1]
                ymin=minmax[2]
                ymax=minmax[3]
                zmin=minmax[4]
                zmax=minmax[5]
                tmin=minmax[6]
                tmax=minmax[7]
                inputvars.append(swim_io.read_nc(copy_from_file,var=v).data[tmin:tmax,zmin:zmax,ymin:ymax,xmin:xmax])
            else: #for now this defaults to all data if minmax=[] (or otherwise)
                inputvars.append(swim_io.read_nc(copy_from_file,var=v).data)
    else:
        inputvars=[]
            
    
    for f in files:
        if verbose:
            print(f)
        add_vars_to_file(f,inputvars,vars2copy,output_prefix)

if __name__ == '__main__':
    filesearch="daily*.nc"
    copy_from_file="/Volumes/G-SAFE/usbr/wrf4km_daily_precip/NARR_04km_OCT2000-SEP2001.nc"
    vars2copy=["XLAT","XLONG"]
    ranges=[[75,150],[95,185]]
    # copy_from_file="/Volumes/G-SAFE/usbr/new_stats/obs/uw.0625/pr/pr.2000.nc"
    # copy_from_file="/Volumes/G-SAFE/usbr/new_stats/obs/maurer.125/pr/21c/nldas_met_update.obs.daily.pr.2000.nc"
    # vars2copy=["lat","lon"]
    # xmin=200;xmax=400
    # xmin=95;xmax=185

    main(filesearch=filesearch,copy_from_file=copy_from_file,ranges=ranges,vars2copy=vars2copy)
    # try:
    #     parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ')
    #     parser.add_argument('filename',action='store')
    #     parser.add_argument('-v', '--version',action='version',
    #             version='Template Parser 1.0')
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
