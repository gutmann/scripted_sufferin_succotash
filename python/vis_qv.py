#!/usr/bin/env python
from __future__ import print_function
import sys
from glob import glob
import datetime
import date_fun

import numpy as np
import matplotlib.pyplot as plt
from stat_down import myio as io

def load_data(filesearch):
    """docstring for load_data"""
    data=io.read_files(filesearch,"qv")
    return data
    
def plt_data(data,files,dates):
    minval=0.0005
    period=365.25
    offset=20
    max_range=0.0055
    max_mean=0.0095
    
    for i in range(len(files)):
        plt.clf();
        data=io.read_nc(files[i],"qv").data
        doy=date_fun.datetime2mjd([dates[i]])-date_fun.date2mjd(dates[i].year,1,1,0,0)
        maxval=np.cos((doy-offset)/period)*max_range+max_mean
        plt.imshow(data[:,0,:],vmin=minval,vmax=maxval,cmap=plt.cm.Blues)
        # else:
        #     plt.imshow(data[i][:,0,:],vmin=minval,vmax=maxval+i*maxvalslope,cmap=plt.cm.Blues)
        plt.colorbar();
        plt.title(dates[i])
        plt.draw()
        plt.savefig("movie/qv_movie_{0:05}.png".format(i))
        # if i%100==0:
            # print(np.round(i/float(len(files))*1000)/10.0,end="% ")
            # print("{0} min={1:8.5f} max={2:8.5f}".format(dates[i],data[:,0,:].min(),data[:,0,:].max()))
            # sys.stdout.flush()
    print("done")

def main():
    """docstring for main"""
    files=glob("output/swim*")
    files.sort()
    startdate=datetime.datetime(2006,1,1,0,0)
    dates=[startdate+datetime.timedelta(i/24.0) for i in range(len(files))]
    plt_data(data=None,files=files,dates=dates)

if __name__ == '__main__':
    main()
    
