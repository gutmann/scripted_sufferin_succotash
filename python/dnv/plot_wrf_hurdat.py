#!/usr/bin/env python

"""
SYNOPSIS

    plot_wrf_hurdat.py [-h] [--verbose] [-v, --version] <filename>

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
import datetime

import matplotlib.pyplot as plt
from mpl_toolkits import basemap

import numpy as np

import mygis

from dnv import load_hurdat, load_wrf_tracks
global verbose
verbose=False

global startdates
startdates = {
    2002:datetime.datetime(2002,10,1,0,0,0),
    2003:datetime.datetime(2003, 7,1,0,0,0),
    2004:datetime.datetime(2004, 8,1,0,0,0),
    2005:datetime.datetime(2005, 7,1,0,0,0),
    2007:datetime.datetime(2007, 9,1,0,0,0),
    # 2008:datetime.datetime(2008, 8,1,0,0,0)}
    2008:datetime.datetime(2008, 7,1,0,0,0)}

global enddates
enddates = {
    2002:datetime.datetime(2002,10,31,23,00,00),
    2003:datetime.datetime(2003, 9,30,23,00,00),
    2004:datetime.datetime(2004, 9,30,23,00,00),
    2005:datetime.datetime(2005,10,31,23,00,00),
    2007:datetime.datetime(2007, 9,30,23,00,00),
    # 2008:datetime.datetime(2008, 8,23,23,00,00)}
    2008:datetime.datetime(2008, 9,30,23,00,00)}

wrf_geo_file="/glade/p/work/gutmann/conus_runs/wrfout_conus_constants.nc"

def plot_wrf_track(d, keys, m, init_date):
    
    geo=mygis.read_geo(wrf_geo_file)
    
    ID  = np.where(keys=="ID")[0]
    lat = np.where(keys=="y")[0]
    lon = np.where(keys=="x")[0]
    time= np.where(keys=="timestep")[0]
    wind= np.where(keys=="wmax")[0]
    size= np.where(keys=="radius")[0]

    if verbose:print(keys)
    if verbose:print(ID, lat, lon, time, wind)
        
    for track in np.unique(d[:,ID]):
        curtrack=np.where(d[:,ID]==track)
        if verbose:print(track)
        if len(curtrack[0])>5:
            if verbose:print("plotting ID: "+str(track))
            x, y = d[curtrack[0],lon].astype(int), d[curtrack[0],lat].astype(int)
            lats = geo.lat[0,y,x+600]
            lons = geo.lon[0,y,x+600]
            
            m.plot(lons,lats,color="black")
            date = init_date + datetime.timedelta(d[curtrack[0][0],time][0]/24.0)
            plt.text(lons[0],lats[0],date.strftime("%Y-%m-%d"))
            
            colors = d[curtrack[0],wind]
            sizes  = d[curtrack[0],size]
            plt.scatter(lons[::6],lats[::6],c=colors[::6],vmin=25,vmax=50, 
                        cmap=plt.cm.jet, s=sizes[::6])

def plot_base_figure(m):
    plt.clf()
    m.bluemarble(alpha=0.5)
    m.drawcoastlines(color="grey")
    m.drawmeridians(np.arange(0,-100,-5), dashes=[5,5], labels=[False, False, False,  True])
    m.drawparallels(np.arange(0,  90, 5), dashes=[5,5], labels=[ True, False, False, False])

def plot_hur_track(data, m, startdate, enddate):
    for i in data:
        if i.strength>2 and (i.dates[0]>=startdate and i.dates[-1]<=enddate):
            m.plot(i.lons,i.lats,"--",color="lightgrey")
            m.scatter(i.lons,i.lats,c=i.winds, vmin=25, vmax=50,cmap=plt.cm.jet, s=i.r64*10)


def main(wrf_filename, outputfile="trackmap_{}.png"):

    if verbose:print("Loading data")
    filename = "hurdat2-2000-2014-060415.txt"
    data = load_hurdat.load(filename)
    keys, wrf_data = load_wrf_tracks.load_data(wrf_filename)

    testyear = int(wrf_filename.split("_")[-1][:4])
    if verbose:print(testyear)
    startdate = startdates[testyear]
    enddate   = enddates[testyear]
    
    
    if verbose:print("Building Map")
    m=basemap.Basemap(llcrnrlon=-100, llcrnrlat=18, urcrnrlon=-69, urcrnrlat=41, resolution="i")

    if verbose:print("Plotting WRF only")
    plot_base_figure(m)
    plot_wrf_track(wrf_data, keys, m, startdate)
    plt.savefig(outputfile.format("wrf_{}".format(testyear)))
    
    if verbose:print("Plotting HURDAT only")
    plot_base_figure(m)
    plot_hur_track(data, m, startdate, enddate)
    plt.savefig(outputfile.format("hurdat_{}".format(testyear)))
    
    
    if verbose:print("Plotting all tracks")
    # plot_base_figure(m)
    plot_wrf_track(wrf_data, keys, m, startdate)
    # plot_hur_track(data, m, startdate, enddate)
    plt.savefig(outputfile.format(testyear))


if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename',nargs="?", action='store', default="ctrl/tracks_subset_2004.txt",
                            help="WRF tracks file")
        parser.add_argument('-o', dest="output_file", action='store', default="trackmap_{}.png",
                            help="Output image filename {} will be replaced by year.")
        parser.add_argument('-v', '--version',action='version',
                version='Template Parser 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose = args.verbose

        exit_code = main(args.filename, args.output_file)
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
