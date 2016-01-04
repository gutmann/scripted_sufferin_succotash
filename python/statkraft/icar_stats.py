#!/usr/bin/env python

"""
SYNOPSIS

    icar_stats.py [-h] [--verbose] [-v, --version] <dir_name>

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

    1.0
"""
from __future__ import absolute_import, print_function, division

import sys
import os
import traceback
import argparse
import glob

import numpy as np
import mygis
from bunch import Bunch
from statkraft import stats

global verbose, start_date, file_search
verbose=False
start_date=0
file_search="icar_????_daily_rain.nc"

def main (dir_name):

    #
    if verbose:
        print('Working on: '+dir_name)
        files=glob.glob(dir_name+file_search)
        files.sort()
        print(files[0]+"  "+files[-1])

    data=mygis.read_files(dir_name+file_search,"precipitation_amount",axis=0)

    dates=np.arange(data.shape[0]+start_date)
    doy=dates%365

    mygis.write(dir_name+"wetfrac.nc",stats.wetfrac(data))
    mygis.write(dir_name+"p99.nc",stats.percentile(data, percentile=0.99))
    mygis.write(dir_name+"mean_annual_precip.nc",stats.mean(data))



if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Generate baseline statistics from icar2daily.py ouput data. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('dir_name',action='store', default="directory to search for data")
        parser.add_argument('-v', '--version',action='version',
                version='ICAR_stats 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose=args.verbose

        exit_code = main(args.dir_name)
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
