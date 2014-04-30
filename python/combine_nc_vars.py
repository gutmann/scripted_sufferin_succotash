#!/usr/bin/env python

"""
SYNOPSIS

    combine_nc_vars.py [-h] [--verbose] [-v, --version] <filename>

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
import mygis as swim_io
# import nio
from netCDF4 import Dataset
verbose=True

def myaverage(data):
    """average all elements of a list together"""
    outputdata=data[0]
    for d in data[1:]:
        outputdata+=d
    return outputdata/float(len(data))

def add_vars_to_file(f1,f2,inputvars,inputatts,varnames,outputdir,output_sub,func,output_prefix):
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
    outputfilename=output_prefix+f1.split("/")[-1].replace(output_sub[0][0],output_sub[0][1])
    # outputfilename=outputfilename.strip(".nc") #strip off the .nc term first... 
    # only strip .nc for nio, netCDF4 expects you to insert the nc
    
    # if there wasn't a directory specified, make sure we aren't specifying "/"
    # if not outputdir:
    #     if len(f1.split("/"))>1:
    #         outputdir="/".join(f1.split("/")[:-1])+"/"
    #     else:
    #         outputdir=""
    # if there wasn't a directory specified, make sure we aren't specifying "/"
    if len(outputdir)==1:
        outputdir=""
    outputf=Dataset(outputdir+outputfilename, 'w',format="NETCDF4")
    inputf1=Dataset(f1, 'r')
    inputf2=Dataset(f2, 'r')
    
    maxlen=0
    # loop over the input variables finding the maximum number of 
    # dimensions we need to know about
    for d in inputvars:
        maxlen=max(maxlen,len(d.shape))
    
    # create an array to hold the references to the dimensions corresponding
    #  to the dimensions of each input variable
    outputdims=np.empty((len(inputvars),maxlen),dtype="S16")
    if len(outputdims.shape)==1:
        outputdims=outputdims.reshape((outputdims.shape[0],1))
    
    dimnum=0
    # loop over the dimensions in the input file
    for i,dim in enumerate(inputf1.dimensions.keys()):
        dimlen=len(inputf1.dimensions[dim]) #in netCDF4python parlence we need a len() around dimension (?)
                                            # for Nio change len above to int
        
        # create that dimension in the outputfile too
        outputf.createDimension(dim,dimlen) #Nio=create_dimension
        # find any input variables that match this dimension
        for j,thisvar in enumerate(inputvars):
            finddim=np.where(np.array(thisvar.shape)==dimlen)
            if len(finddim[0])>0:
                # if we found a match
                # save the dimension in the output dims array
                for k in finddim[0]:
                    outputdims[j,k]=dim
    
    i=0
    # copy the requested variables into the file
    for v,atts,name in zip(inputvars,inputatts,varnames):
        # create the variable
        outputf.createVariable(name,"f",tuple(outputdims[i,:len(v.shape)])) #Nio=create_variable
        # store the data
        outputf.variables[name][:]=v.astype("f")
        outputf.variables[name].setncatts(atts)
        i+=1
        
    # finally, combine the requested variables from the input files to the output file
    v1=output_sub[0][0]
    v2=output_sub[1][0]
    vout=output_sub[0][1]
    thisvar=inputf1.variables[v1]
    var1=thisvar[:]
    var2=inputf2.variables[v2][:]
    
    outputvar=func([var1,var2])
    # create the variable, dimensions and typecode are identical to the input file
    outputf.createVariable(str(vout),thisvar.dtype,thisvar.dimensions) #Nio=create_variable, dtype=typecode()
    # copy the data over
    outputf.variables[str(vout)][:]=outputvar
            
    
    inputf1.close()
    inputf2.close()
    outputf.close()
    
def main (filesearch="*.nc",outputdir="./",output_sub=None,func=None,vars2copy=["lat","lon"],
            ranges=[],output_prefix="",inputvars=None):
    
    files1=glob.glob(filesearch[0])
    files1.sort()
    files2=glob.glob(filesearch[1])
    files2.sort()
    for v in vars2copy:
        ranges.append([])
    
    if func==None:
        func=myaverage
    for f1,f2 in zip(files1,files2):
        if verbose:
            print(f1,f2)
        copy_from_file=f1
        
        inputvars=[]
        inputatts=[]
        for v,minmax in zip(vars2copy,ranges):
            d1=Dataset(copy_from_file)
            print(copy_from_file)
            print(v)
            curatts=dict()
            for k in d1.variables[v].ncattrs():
                curatts[str(k)]=d1.variables[v].getncattr(k)
            inputatts.append(curatts)
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
            
            add_vars_to_file(f1,f2,inputvars,inputatts,vars2copy,outputdir,output_sub,func,output_prefix)
    

if __name__ == '__main__':
    filesearch=["tasmin/*.nc","tasmax/*.nc"]
    outputdir="tas/"
    output_sub=[["tasmin","tas"],["tasmax","tas"]]
    func=myaverage
    
    vars2copy=["lat","lon","time"]

    main(filesearch=filesearch,outputdir=outputdir,output_sub=output_sub,func=func,vars2copy=vars2copy)
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
