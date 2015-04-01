#!/usr/bin/env python

"""
SYNOPSIS

    change_map.py [-h] [-var VARNAME] [-scale SCALE] [-current CURRENT]
                     [-future FUTURE] [-nyears NYEARS] [-geofile GEOFILE] [-v]
                     [--verbose]
                     [filesearch]


DESCRIPTION

    Generate a map of CMIP5 change signal.

    positional arguments:
      filesearch        glob search string for data files (must be escaped)
                        (default: *.nc)

    optional arguments:
      -h, --help        show this help message and exit
      -var VARNAME      name of variable in netcdf files (default: pr)
      -scale SCALE      multiplier scale to be applied to data mean (default:
                        31536000.0)
      -current CURRENT  start date of current period (default: 19900101)
      -future FUTURE    start date of future period (default: 20400101)
      -nyears NYEARS    number of years to compute mean over (default: 10)
      -geofile GEOFILE  file to read geographic data to set up the output map.
                        (default: domain.nc)
      -v, --version     show program's version number and exit
      --verbose         verbose output (default: False)

EXAMPLES

    change_map.py  probabilistic_precipprec_thresh -var data -scale 365 --verbose

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

import matplotlib.pyplot as plt
import numpy as np

import mygis
from stat_down import map_vis

global verbose
verbose=False
global base_year
base_year=1950

def main (filesearch, varname, scale, current, future, nyears, geofile):

    # TODO: Do something more interesting here...
    if verbose:
        print("Verbose output turned on")
        print("Arguments: ")
        print(filesearch, varname, scale, current, future, nyears, geofile)
    
    try:
        if verbose: print("Reading time data")
        times=mygis.read_files(filesearch, "time",axis=0)
    except:
        pass
    
    if verbose: print("Reading primary data")
    data=mygis.read_files(filesearch,varname,axis=0)
    data[data<0]=0
    
    if verbose: print("Reading geographic data")
    geo=mygis.read_geo("domain.nc")
    
    current_year = int(current[:4])
    future_year = int(future[:4])
    
    if verbose: print("Computing current and future means")
    current = scale * data[(current_year-base_year)*365:(current_year-base_year+nyears)*365].mean(axis=0) #*0.9
    future  = scale * data[( future_year-base_year)*365:( future_year-base_year+nyears)*365].mean(axis=0)
    
    
    if verbose: print("Reading map geographic data")
    mapgeo=mygis.read_geo(geofile)
    geobounds=[mapgeo.lat.min(), mapgeo.lat.max(), mapgeo.lon.min(), mapgeo.lon.max()]
    
    if verbose:
        print(future.min(), future.mean(), future.max())
        print(current.min(),current.mean(),current.max())
        delta=future-current
        print(delta.min(),delta.mean(),delta.max())
    
    if verbose: print("Plotting map")
    map_vis.vis((future - current), geo=geobounds, lat=geo.lat, lon=geo.lon, reproject=True,
                cmap=plt.cm.seismic_r, clim=(-300,300), latstep=2, lonstep=5)
    
    plt.savefig("Change_Map.png")
    
    
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Generate a map of CMIP5 change signal. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filesearch',nargs="?",      action='store', default="*.nc",
                            help="glob search string for data files  (must be escaped)")
        parser.add_argument('-var',      dest="varname", action='store', default="pr",
                            help="name of variable in netcdf files")
        parser.add_argument('-scale',    dest="scale",   action='store', default=86400.0*365,type=float,
                            help="multiplier scale to be applied to data mean")
        parser.add_argument('-current',  dest="current", action='store', default="19900101",
                            help="start date of current period")
        parser.add_argument('-future',   dest="future",  action='store', default="20400101",
                            help="start date of future period")
        parser.add_argument('-nyears',   dest="nyears",  action='store', default=10, type=int,
                            help="number of years to compute mean over")
        parser.add_argument('-geofile',  dest="geofile", action='store', default="domain.nc",
                            help="file to read geographic data to set up the output map.")
        
        parser.add_argument('-v', '--version',action='version',
                version='change_map v1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
        
        verbose = args.verbose

        exit_code = main(args.filesearch, args.varname, args.scale, args.current, args.future, args.nyears,args.geofile)
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
