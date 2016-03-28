#!/usr/bin/env python

"""
SYNOPSIS

    template_argparse.py [-h] [--verbose] [-v, --version] <filename>

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

import numpy as np
import matplotlib.pyplot as plt

import mygis
from bunch import Bunch

global verbose
verbose=False

TC_p_threshold =  -2700 # Pa (below max time series p)
land_threshold =    000 # Pa initially set to 101000, but no reason not to track over land too...
wind_threshold =     25 # m/s
dx = 4.0 # km

# permit tracks to move by up to 100 grid cells between tracks
window=100
ntracks=0

def update_tracks(tracks, pd, wind, time):
    # pd is the pressure depression from the timeseries pmax
    newtracks=[]
    
    ny=wind.shape[0]
    nx=wind.shape[1]
    
    for t in tracks[-1]:
        if verbose:print("updating track:{}".format(t.ID))
        ymin=max(t.y-window/2,0)
        ymax=min(t.y+window/2,ny)
        xmin=max(t.x-window/2,0)
        xmax=min(t.x+window/2,nx)
        
        newy,newx = np.unravel_index(np.argmin(pd[ymin:ymax,xmin:xmax]), (ymax-ymin,xmax-xmin))
        wind_window=wind[ymin:ymax,xmin:xmax]
        wmax=wind_window.max()
        newy+=ymin
        newx+=xmin
        pmin=pd[newy,newx]
        # if verbose:print(pmin, TC_p_threshold, wmax, wind_threshold)
        if ((pmin<(TC_p_threshold+1000)) and (wmax>(wind_threshold-10))):
            newtracks.append(get_track_stats(newx,newy, pmin, wmax, wind_window, t.ID, time, xmin, ymin))
            if verbose:print(newtracks[-1])
            if verbose:print("UPDATED!!!", t.ID)
            ymin=max(t.y-window*1.5,0)
            ymax=min(t.y+window*1.5,ny)
            xmin=max(t.x-window*1.5,0)
            xmax=min(t.x+window*1.5,nx)
            pd[ymin:ymax,xmin:xmax]=0
    
    tracks.append(newtracks)
    


def get_track_stats(x,y,pmin,wmax,wind, trackID, time, xmin, ymin):
    
    # find the average radius of winds > wind_threshold
    ny=wind.shape[0]
    npoints=np.zeros((ny,2,2))
    # for each row,
    for i in range(ny):
        # find the first and last x and y locations that exceed the threshold
        tmp=np.where(wind[i]>=wind_threshold)
        if len(tmp[0])>1:
            npoints[i,0,0]=i
            npoints[i,0,1]=tmp[0][0]
            npoints[i,1,0]=i
            npoints[i,1,1]=tmp[0][-1]
    # mask out the rows that did not have any points
    npoints=np.ma.array(npoints,mask=(npoints==0))
    # find the mean center location
    ym=npoints[:,:,0].mean()
    xm=npoints[:,:,1].mean()
    # then compute the distance to each point on the perimeter
    dists=np.sqrt((npoints[:,:,1]-xm)**2+(npoints[:,:,0]-ym)**2)
    # the radius is the mean of those distances
    radius = dists.mean()
    
    # find all points that exceed the threshold so we can compute the area and the mean winds
    high_winds=np.where(wind>wind_threshold)
    
    return Bunch(x=x,y=y, pmin=pmin, wmax=wmax, 
                 radius=radius*dx,
                 area=len(high_winds[0])*(dx**2),
                 wmean=wind[high_winds].mean(), 
                 w3mean=(wind[high_winds]**3).mean(), 
                 ID=trackID, timestep=time,
                 center_x=xm+xmin, center_y=ym+ymin)
    
def find_new_tracks(tracks, pd, wind, time):
    pmin=pd.min()
    if pmin>TC_p_threshold:return # there are no tracks to find
    
    ny=pd.shape[0]
    nx=pd.shape[1]
    newtracks=[]
    global ntracks
    
    while pmin<TC_p_threshold:
        y, x = np.unravel_index(np.argmin(pd),pd.shape)
        ymin=max(y-window/2,0)
        ymax=min(y+window/2,ny)
        xmin=max(x-window/2,0)
        xmax=min(x+window/2,nx)
        
        wmax=wind[ymin:ymax,xmin:xmax].max()
        
        if wmax > wind_threshold:
            if verbose:print("Found new track:{}, {}".format(x,y))
            newtracks.append( get_track_stats(x,y,pmin,wmax,wind[ymin:ymax,xmin:xmax], ntracks, time, xmin, ymin) )
            if verbose:print(newtracks[-1])
            ntracks+=1
            if verbose:print(ntracks)
        
            ymin=max(y-window,0)
            ymax=min(y+window,ny)
            xmin=max(x-window,0)
            xmax=min(x+window,nx)
        pd[ymin:ymax,xmin:xmax]=0
        pmin=pd.min()
    
    tracks[-1].extend(newtracks)
    
    
def write_tracks(tracksequence, filename):
    with open(filename,"w") as f:
        for tracks in tracksequence:
            for t in tracks:
                f.write(str(t)+"\n")
        
def visualize(tracks, pressure, wind, time_step, filename):
    plt.clf()
    
    ax=plt.subplot(1,2,1)
    plt.imshow(wind)
    plt.clim(0,50)
    plt.colorbar()
    plt.title("Wind Speed")
    try:
        for t in tracks:
            plt.text(t.x,t.y, 'x',color="red", horizontalalignment='center', verticalalignment='center')
            if np.isfinite(t.center_x) and np.isfinite(t.center_x):
                circ=plt.Circle((t.center_x,t.center_y), radius=t.radius/dx, color='green', linewidth=2, fill=False)
                ax.add_patch(circ)
                plt.text(t.center_x,t.center_y, 'x',color="green", horizontalalignment='center', verticalalignment='center')
            else:
                circ=plt.Circle((t.x,t.y), radius=t.radius/dx, color='black', linewidth=2, fill=False)
                ax.add_patch(circ)
    except:
        pass
    
    ax=plt.subplot(1,2,2)
    plt.imshow(pressure)
    plt.clim(-5000,0)
    plt.colorbar()
    plt.title("Pressure Anomoly")
    try:
        for t in tracks:
            circ=plt.Circle((t.x,t.y), radius=t.radius/dx, color='black', linewidth=2, fill=False)
            ax.add_patch(circ)
            plt.text(t.x,t.y, 'x',color="red", horizontalalignment='center', verticalalignment='center')
            if np.isfinite(t.center_x) and np.isfinite(t.center_x):
                plt.text(t.center_x,t.center_y, 'x',color="green", horizontalalignment='center', verticalalignment='center')
    except:
        pass
    
    plt.tight_layout()
    plt.savefig(filename+"{:05}.png".format(time_step))
    

def main (filename, outputfilename, p=None, wspd=None, vis_data=False, img_name=""):

    start=0
    end=None
    
    # start=1060
    # end=1100
    if vis_data:
        plt.figure(figsize=(15,6))

    if p==None:
        if verbose: print("Reading PSFC")
        p = mygis.read_nc(filename,"PSFC", returnNCvar=True).data[start:end]
        
    if wspd==None:
        if verbose: print("Reading U10")
        u = mygis.read_nc(filename,"U10", returnNCvar=True).data[start:end]
        if verbose: print("Reading V10")
        v = mygis.read_nc(filename,"V10", returnNCvar=True).data[start:end]
        if verbose:print("Loading wind")
        wspd=np.sqrt(u[:,:600,600:]**2 + v[:,:600,600:]**2)
        del u
        del v
        import gc
        gc.collect()
        
    # subset data to a smaller working region
    # wspd = wspd[:,:600,600:]
    if verbose:print("Loading pressure")
    p=p[:,:600,600:]
    
    if verbose:print("Computing initial values")
    pmax = p.max(axis=0)
    land = (pmax < land_threshold)
    land[:200,400:]=False # set Cuba etc to "ocean" to avoid breaking tracks
    
    tracks=[]
    ntimes=p.shape[0]
    
    for i in range(ntimes):
        if verbose:print("TIME = "+str(i))
        curp = np.ma.array(p[i]-pmax,mask=land)
        if vis_data:
            visp=np.copy(curp)
        curw = np.ma.array(wspd[i],mask=land)
        if i>0:
            update_tracks(tracks, curp, curw, i)
        else:
            tracks.append([])
        
        find_new_tracks(tracks, curp, curw, i)
        if verbose:print("")
        if vis_data:visualize(tracks[-1], visp, curw, i, img_name)
    
    write_tracks(tracks, outputfilename)
        

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Script to track tropical cyclones in the 4km-CONUS WRF output',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename',nargs="?", action='store', default="name of wrf file")
        parser.add_argument('-o',dest="output", action='store', default="tracks.txt")
        parser.add_argument('-f',dest="img_file", action='store', default="map_")
        parser.add_argument ('-p','--plot', action='store_true',
                default=False, help='flag to turn on plot output', dest='plot_data')
        parser.add_argument('-v', '--version',action='version',
                version='Tracker 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose = args.verbose

        exit_code = main(args.filename, args.output, vis_data=args.plot_data, img_name=args.img_file)
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
