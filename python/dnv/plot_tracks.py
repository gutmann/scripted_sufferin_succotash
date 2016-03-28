#!/usr/bin/env python

"""
SYNOPSIS

    plot_tracks.py [-h] [--verbose] [-v, --version] <filename>

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

global verbose
verbose=False

def load_data(filename):
    output=[]
    with open(filename,"r") as f:
        for l in f:
            if len(output)<1:
                if verbose:print("Loading keys")
                keys=[]
                for d in l.split(","):
                    key=d.split(":")[0].replace("'","").strip().replace("{","")
                    keys.append(key)
            
            data=l.split(",")
            values=[]
            for d in data:
                # print(d.split(":")[1])
                try:
                    value=float(d.split(":")[1].replace("'","").strip().replace("{","").replace("}",""))
                except:
                    value=-999
                values.append(value)
                    
            output.append(np.array(values))
    return np.array(keys), np.array(output)
                
def plot_track(d, keys):
    
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
        print(track)
        if len(curtrack[0])>5:
            if verbose:print("plotting ID: "+str(track))
            plt.plot(d[curtrack[0],lon],d[curtrack[0],lat],color="black")
            plt.text(d[curtrack[0],lon][0],d[curtrack[0][0],lat],str(d[curtrack[0][0],time][0]))
            
            colors=d[curtrack[0],wind]
            print(colors.min(), colors.max())
            colors=(colors-25)/30
            colors[colors<0]=0
            colors[colors>0.999]=0.999
            plt.scatter(d[curtrack[0],lon],d[curtrack[0],lat],c=colors, cmap=plt.cm.jet, s=d[curtrack[0],size]*5)
    # if d.shape[0]>5:
    #     for i in range(d.shape[0]):
    #         thistrack=d[i,ID]
    #         if len(np.where(d[:,ID]==thistrack)[0])>5:
    #             c=min(max(0,(d[i,wind]-25)/30),0.99)
    #             # print("")
    #             # print("color=", str(plt.cm.jet(c)[0]))
    #             plt.plot(d[i,lon],d[i,lat],'o',color=plt.cm.jet(c))
    plt.xlim(0,800)
    plt.ylim(0,600)
                 

def main (filename):

    if verbose:print("Loading data")
    keys,data=load_data(filename)
    if verbose:print("Plotting")
    plot_track(data, keys)
    plt.title(filename)
    plt.savefig(filename+".png")
    
if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename',nargs="?", action='store', default="some_file_name")
        parser.add_argument('-f2', dest="file_two", action='store', default="a_Second_filename")
        parser.add_argument('-v', '--version',action='version',
                version='Template Parser 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose = args.verbose

        exit_code = main(args.filename)
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
