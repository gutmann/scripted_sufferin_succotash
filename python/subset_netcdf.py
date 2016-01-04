#!/usr/bin/env python

"""
SYNOPSIS

    subset_netcdf.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION
    Attempt to subset a netcdf file based on the grid coordinates
    Assumes the last coordinate is X and the second to last is Y
    If a variable does not cover the x/y range specified it is output
    without subsetting. 
    
    All variables and attributes are copied over. 

EXAMPLES

    subset_netcdf.py high_res_wrf_output.nc -xmin 320 -xmax 430 -ymin 20 -ymax 150 -o baseline.nc --verbose

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

import mygis
import numpy

from bunch import Bunch

EW_STAGGER_VARS=["U","XLAT_U","XLONG_U","MAPFAC_U","MAPFAC_UX","MAPFAC_UY"]
NS_STAGGER_VARS=["V","XLAT_V","XLONG_V","MAPFAC_V","MAPFAC_VX","MAPFAC_VY","MF_VX_INV"]

def subset_data(ncdata, xmin,xmax, ymin,ymax):
    """docstring for subset_data"""
    sz=ncdata.shape
    if verbose:print("           ",ncdata.shape, xmin,xmax, ymin,ymax)
    if len(sz)<2:return ncdata[:]
    
    if (ymax!=None):
        if (sz[-2]<ymax):
            return ncdata[:]
    else:
        if (sz[-2]<=ymin):
            return ncdata[:]

    if (xmax!=None):
        if (sz[-1]<xmax):
            return ncdata[:]
    else:
        if (sz[-1]<=xmin):
            return ncdata[:]
    
    return ncdata[...,ymin:ymax,xmin:xmax]
        


def main (filename, xmin,xmax, ymin,ymax, outputfile, varnames):
    
    d=mygis.Dataset(filename)
    
    if varnames==None:
        varnames=d.variables.keys()
    else:
        varnames=varnames.split(",")
    
    outputvariables=[]
    if verbose:print("Reading data")
    for v in varnames:
        if verbose:print(v)
        ncdata=d.variables[v]
        
        if (xmax!=None) and (v in EW_STAGGER_VARS):
            xmax+=1
        if (ymax!=None) and (v in NS_STAGGER_VARS):
            ymax+=1
        
        data=subset_data(ncdata, xmin,xmax, ymin,ymax)

        if (xmax!=None) and (v in EW_STAGGER_VARS):
            xmax-=1
        if (ymax!=None) and (v in NS_STAGGER_VARS):
            ymax-=1
        
        atts=Bunch()
        attrlist=ncdata.ncattrs()
        for a in attrlist:
            atts[a]=ncdata.getncattr(a)
        
        outputvariables.append(Bunch(data=data, name=v, attributes=atts, dims=ncdata.dimensions, dtype=data.dtype))
    
    global_atts=Bunch()
    attrlist=d.ncattrs()
    for a in attrlist:
        if a!="history":
            global_atts[a]=d.getncattr(a)
    if "history" in attrlist:
        history="subset by subset_netcdf.py; "+d.getncattr("history")
    else:
        history="subset by subset_netcdf.py"
    
    v=outputvariables[0]
    
    if len(outputvariables)>1:
        outputvariables=outputvariables[1:]
    else:
        outputvariables=None
    
    if verbose:print("Writing output")   
    mygis.write(outputfile,v.data, dtype=v.dtype, varname=v.name, dims=v.dims, attributes=v.attributes, 
                extravars=outputvariables, global_attributes=global_atts)

    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Subset a netcdf file with grid coordinates',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename', action='store')
        parser.add_argument('-xmin', dest="xmin", action='store', default=0, type=int, help="EW grid point min")
        parser.add_argument('-xmax', dest="xmax", action='store', default=0, type=int, help="EW grid point max")
        parser.add_argument('-ymin', dest="ymin", action='store', default=0, type=int, help="NS grid point min")
        parser.add_argument('-ymax', dest="ymax", action='store', default=0, type=int, help="NS grid point max")
        parser.add_argument('-vars', dest="vars", action='store', default=None, help="optional comma seperated list of variables")
        parser.add_argument('-o',    dest="ofile",action='store', default=None, help="outputfile to create")
        parser.add_argument('-v', '--version',action='version',
                version='NetCDF subsetter 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
        
        global verbose
        verbose=args.verbose
        outputfile=args.ofile if args.ofile!=None else "subset_"+args.filename
        
        xmin=args.xmin
        xmax=args.xmax if args.xmax>xmin else None
        ymin=args.ymin
        ymax=args.ymax if args.ymax>ymin else None

        exit_code = main(args.filename, xmin,xmax, ymin,ymax, outputfile, args.vars)
        
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
