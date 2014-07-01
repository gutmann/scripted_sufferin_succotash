#!/usr/bin/env python
from __future__ import print_function
import glob
import sys
import os
import traceback
import argparse

import numpy as np
import scipy.ndimage.morphology as morph

import mygis

# WARNING: these values are all over ridden by argparse, but are defined here to identify module level global variables
#          but CHANGING THE VALUES HERE WILL HAVE NO EFFECT
# Cloud dilation filters
cloud_dilation_doy=160 # day of year to begin expanding the cloud filter.
cloud_threshold=150 # data value to use as threshold, anything greater than this must be "cloud" (actually 100?)

#  Final Filter parameters
nfilter_passes=5    # number of times to loop over the final spatial filter
filter_threshold=50 # number of days of difference permitted between neighbor gridcells
min_reliable_DSD=20 # minimum number of days to snow disapearance to validate this as a "good" gridcell and use it in filtering

# threshold value of SCA to search for
snow_threshold=0

def load_info(pattern="*.tif"):
    """load filenames and read dates from filenames"""
    files=glob.glob(pattern)
    files.sort()
    dates=[int(f[4:7]) for f in files]
    
    return files,dates

def filter_data(data,dates):
    """filter data to expand cloud borders if desired"""
    output=[]
    for date,d in zip(dates,data):
        if (date>cloud_dilation_doy):
            clouds=np.zeros(d.shape)
            clouds[d>cloud_threshold]=1
            clouds=morph.binary_dilation(clouds,iterations=2)
            d[clouds==1]=255
        output.append(d)
    
    return output

def load_data(files,dates,verbose):
    """load data files and apply preliminary filtering"""
    raw_data=[]
    for f in files:
        if verbose:
            print(f,end=", ")
            sys.stdout.flush()
        raw_data.append(mygis.read_tiff(f).data)
        
    if verbose:print("")
    
    return filter_data(raw_data,dates)

def final_spatial_filter(inputdata,verbose):
    """Apply a spatial filter to DSD data to remove spurious points"""
    
    # window size to use in the filter
    wsize=3
    whalf=(wsize-1)/2
    #copy the input data so we don't alter the array as we work on it
    data=inputdata.copy()
    for i in range(whalf,data.shape[0]-whalf):
        if verbose:
            if (i%int(round((data.shape[0]-whalf*2)/10.0)))==0:
                print((100*i)/(data.shape[0]-whalf*2),end="% ")
                sys.stdout.flush()
                
        for j in range(whalf,data.shape[1]-whalf):
            # set up the local area of interest
            local_window=inputdata[i-whalf:i+whalf,j-whalf:j+whalf]
            # find the points that are reliable
            local_window=local_window[local_window>min_reliable_DSD]
            # if there are reliable points
            if local_window.size>0:
                # find the minimum value of the surrounding reliable points
                local_min=local_window.min()
                # if the current point is greate than this value by more than a given threshold
                # replace the current point with the median of the surrounding good points
                if (data[i,j]-local_min)>filter_threshold:
                    data[i,j]=np.median(local_window[local_window<data[i,j]])
    if verbose:print("")
    return data
            

def main(pattern="*.tif",ndays=2,verbose=True):
    files,dates= load_info(pattern)
    
    if verbose:print("Loading data")
    data=load_data(files,dates,verbose=verbose)
    
    # set up variables we will need
    last_snow=np.zeros(data[0].shape)
    snow_on_again=np.zeros(data[0].shape)
    snow_on=np.empty(data[0].shape,dtype=bool)
    snow_on[:]=True
    
    last_date=0
    for date,img in zip(dates,data):
        if verbose:
            print(date,end=", ")
            sys.stdout.flush()
        snow_off_now=(img==snow_threshold)
        #points to mark as snow off are those that had snow previously, but don't now
        snow_off_points=np.where(snow_on & snow_off_now)
        
        # at those points, set DSD as occuring between the current date and the last snow covered date
        last_snow[snow_off_points]=(date+last_date)/2.0
        # also mark them as no longer having snow
        snow_on[snow_off_points]=False
        # and set the snow on again counter to 0
        snow_on_again[snow_off_points]=0

        # at points that have snow, add one to the snow on again counter
        snow_on_again[(img<cloud_threshold)&(img>snow_threshold)]+=1
        # after more than n days in a row set snow_on flag again
        snow_on[snow_on_again>ndays]=True
        last_date=date
    
    if verbose:print("")
    
    for i in range(nfilter_passes):
        print("Filtering data: pass {}".format(i))
        outputdata=final_spatial_filter(last_snow,verbose)
        last_snow[:]=outputdata[:]
    mygis.write("DSD.nc",outputdata)
    
    

if __name__ == '__main__':
    global nfilter_passes
    global filter_threshold
    global min_reliable_DSD
    global cloud_dilation_doy
    global cloud_threshold
    global snow_threshold

    try:
        parser= argparse.ArgumentParser(description='Find the date of snow disappearance in a series of SCA images. ')

        parser.add_argument('search',nargs="?",action='store',default="*.tif",
                            help="glob string to use when searching for filenames, be sure to properly \escape")
        parser.add_argument('-ndays',type=int,nargs="?", dest='ndays_snow_on',action='store',default=2,
                            help="Number of days in a row with snow on to reset snow actually on")
        parser.add_argument('-nfilters',type=int,nargs="?", dest='nfilter_passes',action='store',default=5,
                            help="Number of spatial filter passes to run")
        parser.add_argument('-filter',type=int,nargs="?", dest='filter_threshold',action='store',default=50,
                            help="Threshold to use for allowable spatial difference in DSD")
        parser.add_argument('-minDSD',type=int,nargs="?", dest='min_reliable_DSD',action='store',default=20,
                            help="Minimum DSD to treat as reliable for use in the filter")
        parser.add_argument('-cloud_threshold',type=int,nargs="?", dest='cloud_threshold',action='store',default=150,
                            help="Threshold to use to define cloud or other bad data")
        parser.add_argument('-cloud_day',type=int,nargs="?", dest='cloud_dilation_doy',action='store',default=160,
                            help="Day of the year to begin dilating bad (cloudy) data")
        parser.add_argument('-snow_threshold',type=int,nargs="?", dest='snow_threshold',action='store',default=0,
                            help="Threshold to use to define snow cover [default=0]")
        parser.add_argument('-v', '--version',action='version',
                            version='calc_dsd v1.0')
        parser.add_argument ('--verbose', action='store_true', default=False, 
                            help='verbose output', dest='verbose')
        args = parser.parse_args()

        nfilter_passes=args.nfilter_passes
        filter_threshold=args.filter_threshold
        min_reliable_DSD=args.min_reliable_DSD
        cloud_dilation_doy=args.cloud_dilation_doy
        cloud_threshold=args.cloud_threshold
        
        snow_threshold=args.snow_threshold

        exit_code = main(args.search, args.ndays_snow_on, args.verbose)
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
    
    