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
import datetime as dt

import numpy as np
import mygis
from bunch import Bunch
from statkraft import stats

global verbose, rain_file_search, tmin_file_search, tmax_file_search, days_per_year
verbose=False
rain_file_search="icar_????_daily_rain.nc"
tmin_file_search="icar_????_daily_tmin.nc"
tmax_file_search="icar_????_daily_tmax.nc"
days_per_year=365

global month_start, month_end, calendar, start_date

start_date=dt.datetime(1979,01,01,00,00)
calendar="gregorian"

month_start = np.array([ 0, 31, 59, 90, 120,151,181,212,243,273,304,334])
month_end   = np.array([31, 59, 90, 120,151,181,212,243,273,304,334,365])


def compute_dates(ndays):
    if calendar=="noleap":
        dates=[dt.datetime(2001,1,1)+dt.timedelta(i) for i in range(ndays)]
    elif calendar=="gregorian":
        dates=[start_date+dt.timedelta(i) for i in range(ndays)]
    return np.array(dates)

def rain_stats(dir_name):
    if verbose:
        files=glob.glob(dir_name+rain_file_search)
        files.sort()
        print(files[0]+"  "+files[-1])

    data=mygis.read_files(dir_name+rain_file_search,"precipitation_amount",axis=0, verbose=verbose)

    if verbose: print("Computing wet day fraction")
    mygis.write(dir_name+"annual_wetfrac.nc",stats.wetfrac(data))

    if verbose: print("99th percentile")
    mygis.write(dir_name+"annual_p99.nc",stats.percentile(data, percentile=0.99))

    if verbose: print("mean annual precip")
    mygis.write(dir_name+"annual_mean.nc",stats.mean(data)*days_per_year)

    ndays=data.shape[0]
    dates=compute_dates(ndays)
    curmonth=np.empty(ndays,dtype=bool)

    for month in range(12):
        for i in range(ndays):
            curmonth[i] = (dates[i].month == (month+1))


        if verbose: print("month{:02}".format(month+1))
        if verbose: print("     wet day fraction")
        mygis.write(dir_name+"month{:02}_wetfrac.nc".format(month+1),stats.wetfrac(data[curmonth,:,:]))

        if verbose: print("     99th percentile")
        mygis.write(dir_name+"month{:02}_p99.nc".format(month+1),stats.percentile(data[curmonth,:,:], percentile=0.99))

        if verbose: print("     mean precip")
        mygis.write(dir_name+"month{:02}_mean_precip.nc".format(month+1),stats.mean(data[curmonth,:,:])*(month_end[month]-month_start[month]))

    return data

def tmax_stats(dir_name):
    if verbose:
        files=glob.glob(dir_name+tmax_file_search)
        files.sort()
        print(files[0]+"  "+files[-1])

    data=mygis.read_files(dir_name+tmax_file_search,"daily_maximum_temperature",axis=0, verbose=verbose)

    if verbose: print("mean annual tmax")
    mygis.write(dir_name+"annual_mean_tmax.nc",stats.mean(data))

    ndays=data.shape[0]
    dates=compute_dates(ndays)
    curmonth=np.empty(ndays,dtype=bool)

    for month in range(12):
        for i in range(ndays):
            curmonth[i] = (dates[i].month == (month+1))

        if verbose: print("month{:02}".format(month+1))
        if verbose: print("     mean")
        mygis.write(dir_name+"month{:02}_mean_tmax.nc".format(month+1),stats.mean(data[curmonth,:,:]))

    return data

def tmin_stats(dir_name):
    if verbose:
        files=glob.glob(dir_name+tmin_file_search)
        files.sort()
        print(files[0]+"  "+files[-1])

    data=mygis.read_files(dir_name+tmin_file_search,"daily_minimum_temperature",axis=0, verbose=verbose)

    if verbose: print("mean annual tmin")
    mygis.write(dir_name+"annual_mean_tmin.nc",stats.mean(data))

    ndays=data.shape[0]
    dates=compute_dates(ndays)
    curmonth=np.empty(ndays,dtype=bool)

    for month in range(12):
        for i in range(ndays):
            curmonth[i] = (dates[i].month == (month+1))

        if verbose: print("month{:02}".format(month+1))
        if verbose: print("     mean")
        mygis.write(dir_name+"month{:02}_mean_tmin.nc".format(month+1),stats.mean(data[curmonth,:,:]))

    return data

def dtr_stats(dir_name, tmax_data,tmin_data):
    data=tmax_data-tmin_data

    if verbose: print("mean annual diurnal temperature range")
    mygis.write(dir_name+"annual_mean_dtr.nc",stats.mean(data))

    ndays=data.shape[0]
    dates=compute_dates(ndays)
    curmonth=np.empty(ndays,dtype=bool)

    for month in range(12):
        for i in range(ndays):
            curmonth[i] = (dates[i].month == (month+1))

        if verbose: print("month{:02}".format(month+1))
        if verbose: print("     mean")
        mygis.write(dir_name+"month{:02}_mean_dtr.nc".format(month+1),stats.mean(data[curmonth,:,:]))




def main (dir_name, run_rain=False, run_tmin=False, run_tmax=False):

    if dir_name[-1]!="/":dir_name+="/"
    if verbose: print('Working on: '+dir_name)

    if run_rain:
        rain=rain_stats(dir_name)

    if run_tmin:
        tmin=tmin_stats(dir_name)

    if run_tmax:
        tmax=tmax_stats(dir_name)

    if run_tmax and run_tmin:
        dtr=dtr_stats(dir_name, tmin, tmax)


if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Generate baseline statistics from icar2daily.py ouput data. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('dir_name',action='store', help="directory to search for data")
        parser.add_argument('-v', '--version',action='version',
                version='ICAR_stats 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        parser.add_argument ('--rain', action='store_true',
                default=False, help='calculate precip stats', dest='rain')
        parser.add_argument ('--tmax', action='store_true',
                default=False, help='calculate tmax stats', dest='tmax')
        parser.add_argument ('--tmin', action='store_true',
                default=False, help='calculate tmin stats', dest='tmin')
        parser.add_argument ('--noleap', action='store_true',
                default=False, help='use a noleap calendar', dest='noleap')
        parser.add_argument ('--gregorian', action='store_true',
                default=False, help='use a gregorian calendar', dest='gregorian')

        args = parser.parse_args()

        if args.noleap:
            calendar="noleap"
        elif args.gregorian:
            calendar="gregorian"
        else:
            if args.dir_name.split("_")[2][:3]=="era":
                calendar="gregorian"
            else:
                calendar="noleap"

        verbose=args.verbose
        if verbose:print("Using a {} calendar".format(calendar))

        exit_code = main(args.dir_name, run_rain=args.rain, run_tmin=args.tmin, run_tmax=args.tmax)
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
