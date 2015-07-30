#!/usr/bin/env python

"""
SYNOPSIS

    subset_netcdf.py [-h] [--verbose] [-v, --version] <filename>

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

import mygis
import numpy

from bunch import Bunch

def subset_data(ncdata, xmin,xmax, ymin,ymax):
    """docstring for subset_data"""
    sz=ncdata.shape
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
    for v in varnames:
        print(v)
        ncdata=d.variables[v]
        data=subset_data(ncdata, xmin,xmax, ymin,ymax)
        
        atts=Bunch()
        attrlist=ncdata.ncattrs()
        for a in attrlist:
            atts[a]=ncdata.getncattr(a)
        
        # e.data,e.name,e.dims,e.dtype,e.attributes
        outputvariables.append(Bunch(data=data, name=v, attributes=atts, dims=ncdata.dimensions, dtype='f'))
    
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
