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
# make pretty CONUS wide WRF movies

import sys
import os
import traceback
import argparse
import datetime
import glob

import mygis
import matplotlib.pyplot as plt
import numpy as np

import custom_cmap
from mpl_toolkits.basemap import Basemap

def map_vis(data,geo=[],title="",vmin=None,vmax=None,cmap=None,showcolorbar=True,
            latlabels=[1,0,0,0],lonlabels=[0,0,0,1],cbar_label=None,m=None,imgdata=None):
    """Plot a map of data using the bounds in geo=[lllat,urlat,lllon,urlon]
    
    Optionally specify a map title, min and max value and colormap
    """
    if m==None:
        raise ValueError("Must input a basemap instance")
    if imgdata==None:
        mapimg=m.imshow(data,vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        mapimg=m.imshow(imgdata)
    
    m.drawparallels(np.arange(20,60,5.),labels=latlabels,dashes=[1,4])
    m.drawmeridians(np.arange(-130,-65,10.),labels=lonlabels,dashes=[1,4])
    m.drawstates(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    
    if showcolorbar:
        cb=m.colorbar()
        if cbar_label!=None:
            cb.set_label(cbar_label)
    
    if title:
        plt.title(title,x=0.05,ha="left")
        
    return mapimg

    
def plt_data(data, curdate, m, overlay, map_img=None, overimg=None,imgdata=None):
    minval=0.0005
    period=365.25
    offset=20
    max_range=0.005
    # max_mean=0.0095
    max_mean=0.015
    
    
    doy=curdate - datetime.datetime(2000,1,1,0,0,0)
    doy=(doy.days+doy.seconds/86400.0)%365.25
    print(str(curdate))
    maxval=np.cos((doy-offset-period/2)/period*2*np.pi)*max_range+max_mean
    
    if map_img==None:
        # print("Slow plot")
        plt.clf();
        map_img=map_vis(data,vmin=minval,vmax=maxval,cmap=plt.cm.Blues,title=str(curdate),
                showcolorbar=imgdata==None, cbar_label="Water Vapor [kg/kg]",m=m, imgdata=imgdata)
        if (overlay!=None) and (imgdata==None):
            overimg=m.imshow(overlay)
        # map_img.set_cmap(plt.cm.Blues)
    else:
        if imgdata==None:
            map_img.set_data(data)
            map_img.set_clim((minval,maxval))
            if overimg!=None:
                overimg.set_data(overlay)
        else:
            # print("Fast plot")
            map_img.set_data(imgdata)
        plt.title(str(curdate),x=0.05,ha="left")
            
    plt.draw()
    
    return map_img,overimg

def main (data_dir, year, month, day, hour):
    
    plt.figure(figsize=(13,10))
    files=glob.glob(data_dir+"/wrf2d*00")
    files.sort()
    print(data_dir)
    print(files[0])
    # files=files[:10]
    # files=files[960:]
    dx=4000
    test_data=mygis.read_nc(files[0],"Q2",returnNCvar=True)
    nx=test_data.data.shape[2]
    ny=test_data.data.shape[1]
    nt=test_data.data.shape[0]
    test_data.ncfile.close()
    
    m = Basemap(width=nx*dx,height=ny*dx,
                rsphere=(6378137.00,6356752.3142),\
                resolution='l',area_thresh=10000.,projection='lcc',\
                lat_1=28.,lat_2=50.,lat_0=39.7,lon_0=-98.0)
    
    start_date=datetime.datetime(year,month,day,hour)
    ntimes=len(files)*nt
    times_per_day=24.0
    # print("Loading water vapor data")
    # data=mygis.read_files(files,"Q2",axis=0)
    # print("Loading precip data")
    # precipdata=mygis.read_files(files,"PREC_ACC_NC",axis=0)
    
    precip_cmap=custom_cmap.subset(plt.cm.jet,clim=(130,255))
    max_precip=30.0
    min_precip=1
    map_img,overimg=None,None
    #qv color scale
    minval=0.0005
    period=365.25
    offset=20
    max_range=0.005
    # max_mean=0.0095
    max_mean=0.015
    
    savetime=0.0
    plottime=0.0
    readtime=0.0
    setuptime=0.0
    import time
    print("Generating movie for date:")
    curtimes=0
    nexttime=0
    imgdata=None
    overlay=None
    full_runt0=time.time()
    for i in range(ntimes):
        curdate=start_date + datetime.timedelta(i/times_per_day)
        if i>=nexttime:
            t0=time.time()
            # print("Reading data")
            fulldata=mygis.read_nc(files[i],"Q2").data
            curtimes=fulldata.shape[0]
            
            fullprecipdata=mygis.read_nc(files[i],"PREC_ACC_NC").data
            readtime+=(time.time()-t0)
            lasttime=i
            nexttime=i+curtimes
            
        
        doy=curdate - datetime.datetime(2000,1,1,0,0,0)
        doy=(doy.days+doy.seconds/86400.0)%365.25
        maxval=np.cos((doy-offset-period/2)/period*2*np.pi)*max_range+max_mean
        
        
        t0=time.time()
        data=fulldata[i-lasttime]
        precipdata=fullprecipdata[i-lasttime]
        
        overlay=precip_cmap(precipdata/max_precip) # create an overlay for precipitation
        # alpha=np.zeros(overlay.shape[:-1])
        # alpha=precipdata/(max_precip*3)+0.8
        # alpha[alpha>1]=1
        precip_mask=precipdata<min_precip
        # alpha[precip_mask]=0
        overlay[:,:,3]=1
        overlay[:,:,3]=precipdata/(max_precip*3)+0.8 # set alpha to 75% everywhere
        overlay[:,:,3][overlay[:,:,3]>1]=1
        overlay[:,:,3][precip_mask]=0 # set alpha to 0 if precip less than x
        
        imgdata = plt.cm.Blues(data/maxval)
        for rgb in range(3): imgdata[:,:,rgb][~precip_mask] = overlay[:,:,rgb][~precip_mask]
        # for rgb in range(3):
        #     imgdata[:,:,rgb] = (overlay[:,:,rgb]*(1-alpha)) + (overlay[:,:,rgb] * alpha)
        setuptime+=time.time()-t0
        
        
        t0=time.time()
        # print("Plotting data")
        # map_img,overimg=None,None
        map_img,overimg = plt_data(data, curdate, m,overlay,map_img,overimg,imgdata=imgdata)
        # overimg.set_cmap(precip_cmap)
        # plt.savefig("movie/qv_movie_{0:05}.png".format(i))
        if i>0:
            plottime+=(time.time()-t0)
        else:
            plot_init=time.time()-t0
        t0=time.time()
        # print("Saving image")
        plt.savefig("movie/qv_movie_{}.png".format(i))#,dpi=50)
        savetime+=(time.time()-t0)
    
    print("Finished:"+str(time.time()-full_runt0))
    print("Reading:"+str(readtime))
    print("Settup:"+str(setuptime))
    print("Ploting:"+str(plottime)+" "+str(plot_init))
    print("Saving:"+str(savetime))
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('data_dir',nargs="?", action='store', default="data")
        parser.add_argument('-Y',  dest="year", action='store', default=2005,type=int, help="start date year")
        parser.add_argument('-M',  dest="month",action='store', default=10,  type=int, help="start date month")
        parser.add_argument('-D',  dest="day",  action='store', default=1,   type=int, help="start date day")
        parser.add_argument('-hr', dest="hour", action='store', default=0,   type=int, help="start date hour")
        parser.add_argument('-v', '--version',action='version',
                version='qv_movie 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
        
        global verbose
        verbose=args.verbose
        
        exit_code = main(args.data_dir, args.year, args.month, args.day, args.hour)
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
