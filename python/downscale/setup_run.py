#!/usr/bin/env python

"""
SYNOPSIS

    usage: setup_run.py [-h] [-n [NPROC]] [-s [STATION]] [-t [TEMPLATE]]
                        [-d [STATION_DIR]] [-v] [--verbose]
    
    Set up config files for a parallel downscaling run.
    
    optional arguments:
      -h, --help        show this help message and exit
      -n [NPROC]        Number of ways to split up the station file (number of
                        processors to use). (default: 32)
      -s [STATION]      Filename of file listing all stations/grid points to
                        downscale (default: stnlist)
      -t [TEMPLATE]     Filename of a file to use as a template when creating new
                        config files (default: config_template)
      -d [STATION_DIR]  Directory containing the station data (and station file)
                        (default: maurer-sites/)
      -v, --version     Print version
      --verbose         verbose output (default: False)
      
      
DESCRIPTION

    Takes a master station list file and divides it evenly into n-station files. 
    Also sets up a downscale config file for each stationlist file based on a template file. 
    
EXAMPLES

    setup_run.py -n 32 -s stnlist -t template_file -d conus_sites/ --verbose

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
import glob

from bunch import Bunch
import numpy as np

def read_station_file(filename):
    """read_station_file example file
    NSITES 3
        #STNID    LAT         LONG   ELV_FT   STN_NAME
    0000000maurer, 34.0625, -113.938, 0000, Grid_point
    0000112maurer, 34.1875, -113.938, 0000, Grid_point
    0000224maurer, 34.3125, -113.938, 0000, Grid_point
    """
    with open(filename,"rU") as f:
        # extract the number of sites from the line "NSITES 3" (for example)
        n=int(f.next().split()[1])
        header=f.next()
        output=[]
        for l in f:
            # there have to be at least 5 columns of data with 4 column separators
            if len(l.strip())>=9:
                output.append(l)
            else:
                print("ERROR: line does not contain a valid station entry. ")
                print(l)
        
    
    if (n!=len(output)):
        print("ERROR number of sites defined in file does not match number of sites lists")
        print("  Using the actual number of sites defined:"+str(len(output)))
        print("  instead of the number defined:"+str(n))
        n=len(output)
    
    return Bunch(n=n, data=output)
        
def split_stations(n, stationlist):
    """evenly split a list of stations into n lists"""
    output=[]
    nsites=int(np.ceil(stationlist.n/n))
    for i in range(n):
        if (i+1*nsites)>stationlist.n:
            output.append(Bunch(data=stationlist.data[i*nsites:],name=str(i+1)))
        else:
            output.append(Bunch(data=stationlist.data[i*nsites:(i+1)*nsites],name=str(i+1)))
    return output

def write_station_file(station, basefile):
    """docstring for write_station_file"""
    outputfile="stn_"+station.name+basefile.split("/")[-1]
    output_dir="/".join(basefile.split("/")[:-1])
    if len(output_dir)>1:
        if output_dir[-1]!="/":output_dir+="/"
    
    with open(output_dir+outputfile,"w") as f:
        f.write("NSITES {}\n".format(len(station.data)))
        f.write("#STNID    LAT         LONG   ELV_FT   STN_NAME\n")
        f.writelines(station.data)
    
    return output_dir+outputfile

def write_config_file(template, station, name,verbose):
    f=open(template,"rU")
    config_file="config_"+name
    if verbose:print(config_file)
    with open(config_file,"w") as o:
        for l in f:
            l=l.replace("__STATION_FILE__",station)
            l=l.replace("__OUTPUT_FILE__","stn_"+name+"_prcp_norm.nc")
            o.write(l)
    f.close()
    return config_file

def main (nways, nsites, station_file, template_file, directory, verbose):

    if verbose:
        print(nways, nsites, station_file, template_file, directory)
    
    if verbose:print("Reading station file: "+directory+station_file)
    stations=read_station_file(directory+station_file)
    
    if (nsites>0):
        nways=int(np.ceil(stations.n/nsites))
    if verbose:print("N splits = "+str(nways))
    
    divided_stations=split_stations(nways, stations)
    if verbose:print("Nstations:"+str(len(divided_stations)))
    
    with open("mpi_appfile","w") as f:
        for s in divided_stations:
            curstn_file=write_station_file(s, directory+station_file)
            if verbose:print(curstn_file)
            config_file=write_config_file(template_file, curstn_file, s.name,verbose)
            f.write("-np 1 downscale {0} \n".format(config_file))
            # f.write("-np 1 downscale {0} &>{0}.out \n".format(config_file))
    
    
if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Set up config files for a parallel downscaling run. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-n',dest='nproc',nargs='?',action='store',default=32, type=int, 
                            help="Number of ways to split up the station file (number of processors to use). ")
        parser.add_argument('-nsites',dest='nsites',nargs='?',action='store',default=-1, type=int, 
                            help="Number of Stations to put into each file, effectively overrides -n")
        parser.add_argument('-s',dest='station',nargs='?',action='store',default="stnlist", 
                            help="Filename of file listing all stations/grid points to downscale")
        parser.add_argument('-t',dest='template',nargs='?',action='store',default="config_template",
                            help="Filename of a file to use as a template when creating new config files")
        parser.add_argument('-d',dest='station_dir',nargs='?',action='store',default="maurer-sites/",
                            help="Directory containing the station data (and station file)")
        parser.add_argument('-v', '--version',action='version',
                            version='1.0', help="Print version")
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        parser.add_argument ('--config', action='store_true',
                default=False, help='Just write config files from the template', dest='config')
        args = parser.parse_args()

        exit_code=None
        if args.config:
            stnfiles=glob.glob(args.station_dir+args.station)
            print(len(stnfiles))
            print(args.station_dir+args.station)
            for s in stnfiles:
                write_config_file(args.template, s, s.split("_")[-1], args.verbose)
        else:
            exit_code = main(args.nproc, args.nsites, args.station, args.template, args.station_dir, args.verbose)
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
