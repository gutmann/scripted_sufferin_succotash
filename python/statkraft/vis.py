#!/usr/bin/env python

"""
SYNOPSIS

    vis.py [-h] [--verbose] [-v, --version] <data_dir>, <icar_dir>

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
from __future__ import absolute_import, print_function, division

import sys
import os
import traceback
import argparse
import glob

import matplotlib.pyplot as plt
import numpy as np

import mygis
from statkraft import map, icar, obs

global verbose
verbose=False

base_icar_file = "/glade/p/work/gutmann/statkraft/icar_output/icar_run_era/icar_2000_daily_tmin.nc"
base_obs_file  = "/glade/p/work/gutmann/statkraft/data/annual/PREC_2000.nc"

file_search="{time}_mean_{var}.nc"
wetfracfile_search="{time}_wetfrac.nc"

cmaps=dict(precip=plt.cm.gist_rainbow, temp=plt.cm.gist_rainbow_r, wetfrac=plt.cm.gist_rainbow)
clims=dict(mprecip=(0,400),aprecip=(0,5000), atemp=(-7,12),mtemp=(-10,17), awetfrac=(0,1),mwetfrac=(0,1))
# clims=dict(precip=None, temp=None)
title="{}  {}  {}"
figure_file="{}_{}.png"

def plot_data(obs,icar,time,var_name, output_dir, m=None):
    plt.figure(figsize=(8,5))

    plt.subplot(1,2,1)
    data=np.ma.masked_array(obs.data,mask=obs.data<-100)
    lat=obs.lat
    lon=obs.lon
    m=map.vis(data, lat, lon, cmap=cmaps[var_name], clim=clims[time[0]+var_name], title=title.format("Obs",time,var_name),m=m)

    plt.subplot(1,2,2)
    data=np.ma.masked_array(icar.data,mask=icar.data<-100)
    lat=icar.lat
    lon=icar.lon
    m=map.vis(data, lat, lon, cmap=cmaps[var_name], clim=clims[time[0]+var_name], title=title.format("ICAR",time,var_name),m=m)

    plt.tight_layout()
    plt.savefig(output_dir+figure_file.format(time,var_name))
    plt.close()
    return m

def plot_wetfrac(time, obs_dir, icar_dir, output_dir, igeo,ogeo,m=None):
    icar_file = glob.glob(icar_dir+wetfracfile_search.format(time=time))[0]
    if verbose:print("Reading:"+icar_file)
    icar_data = mygis.read_nc(icar_file).data

    obs_file = glob.glob(obs_dir+wetfracfile_search.format(time=time))[0]
    if verbose:print("Reading:"+obs_file)
    obs_data=mygis.read_nc(obs_file).data

    igeo.data=icar_data
    ogeo.data=obs_data

    if verbose:print("Plotting wetfrac")
    m=plot_data(ogeo,igeo,time,"wetfrac",output_dir,m=m)


def main (obs_dir,icar_dir, output_dir):

    igeo=icar.geo(base_icar_file)
    ogeo=obs.geo(base_obs_file, subset=True)
    times=["annual"]
    times.extend(["month{:02}".format(month) for month in range(1,13)])
    m=None
    for time in times:
        for variable in ["temp","precip"]:
            if verbose:print(time+"  "+variable)
            if variable=="temp":
                icar_file = glob.glob(icar_dir+file_search.format(time=time,var="tmin"))[0]
                if verbose:print("Reading:"+icar_file)
                icar_tmin = mygis.read_nc(icar_file).data
                icar_file = glob.glob(icar_dir+file_search.format(time=time,var="tmax"))[0]
                if verbose:print("Reading:"+icar_file)
                icar_tmax = mygis.read_nc(icar_file).data
                icar_data = (icar_tmin+icar_tmax)/2-4

            else:
                icar_file = glob.glob(icar_dir+file_search.format(time=time,var=variable))[0]
                if verbose:print("Reading:"+icar_file)
                icar_data = mygis.read_nc(icar_file).data

            obs_file = glob.glob(obs_dir+file_search.format(time=time,var=variable))[0]
            if verbose:print("Reading:"+obs_file)
            obs_data=mygis.read_nc(obs_file).data

            igeo.data=icar_data
            ogeo.data=obs_data

            if verbose:print("Plotting")
            m=plot_data(ogeo,igeo,time,variable,output_dir,m=m)

        plot_wetfrac(time,obs_dir,icar_dir, output_dir,igeo,ogeo,m=m)



if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Visualize ERA forced ICAR output as compared to met.no obs. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('data_dir', nargs="?", action='store', default="data/annual/")
        parser.add_argument('icar_dir', nargs="?", action='store', default="icar_output/icar_run_era/")
        parser.add_argument('-o', nargs="?", action='store', default="./", dest="output")
        parser.add_argument('-v', '--version',action='version',
                version='vis 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose = args.verbose

        exit_code = main(args.data_dir, args.icar_dir, args.output)
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
