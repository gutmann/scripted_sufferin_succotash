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
      
DESCRIPTION

    Takes a master station list file and divides it evenly into n-station files. 
    Also sets up a downscale config file for each stationlist file based on a template file. 
    
EXAMPLES

    setup_run.py -n 32 -s stnlist -t template_file -d conus_sites/ --verbose

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

def main (station_file, template_file, directory, verbose):

    # TODO: Do something more interesting here...
    print(station_file, template_file, directory, verbose)
    
if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Set up config files for a parallel downscaling run. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-n',dest='nproc',nargs='?',action='store',default=32, type=int, 
                            help="Number of ways to split up the station file (number of processors to use). ")
        parser.add_argument('-s',dest='station',nargs='?',action='store',default="stnlist", 
                            help="Filename of file listing all stations/grid points to downscale")
        parser.add_argument('-t',dest='template',nargs='?',action='store',default="config_template",
                            help="Filename of a file to use as a template when creating new config files")
        parser.add_argument('-d',dest='station_dir',nargs='?',action='store',default="maurer-sites/",
                            help="Directory containing the station data (and station file)")
        parser.add_argument('-v', '--version',action='version',
                version='1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        exit_code = main(args.station, args.template, args.station_dir, args.verbose)
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
