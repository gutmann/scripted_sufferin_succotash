#!/usr/bin/env python
from __future__ import print_function
import sys
from glob import glob
import datetime

import numpy as np
import matplotlib.pyplot as plt
from stat_down import myio as io

def load_data(filesearch):
    """docstring for load_data"""
    data=io.read_files(filesearch,"qv")
    return data
    
def plt_data(data,dates):
    minval=0.0005
    maxval=0.004
    maxvalslope=1.1e-6
    for i in range(len(data)):
        plt.clf();
        plt.imshow(data[i][:,0,:],vmin=minval,vmax=maxval+i*maxvalslope,cmap=plt.cm.Blues)
        plt.colorbar();
        plt.title(dates[i])
        plt.draw()
        plt.savefig("movie_simple/qv_movie_{0:05}.png".format(i))
        if i%100==0:
            print(np.round(i/float(len(data))*1000)/10.0,end="% ")
            sys.stdout.flush()
    print("done")

def main():
    """docstring for main"""
    data=load_data("output_simple/swim*")
    startdate=datetime.datetime(2008,1,1,0,0)
    dates=[startdate+datetime.timedelta(i/24.0) for i in range(len(data))]
    plt_data(data,dates)

if __name__ == '__main__':
    main()
    
