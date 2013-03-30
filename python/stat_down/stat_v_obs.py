#!/usr/bin/env python
from __future__ import print_function

"""
SYNOPSIS

    stat_v_obs.py [-h] [--verbose] [-v, --version] [control_file]

DESCRIPTION

    Compares all statistical downscaling methods to observations

    Inputs:
        Optional control_file which specifies all other files to be read/compared.  
            default is the "standard" set of USBR-USACE statistical downscaling producs
            and the UW-observations dataset

    Outputs: 
        For each statistical downscaling method output
        RMS and Bias maps (png and netCDF) and time-series (txt and png plots)

EXAMPLES

    stat_v_obs.py stat_info.txt
    
    short_run.txt example :
        obs = obs/uw.4km/ : *.200[1-2].nc
        datasets = narr, ncep
        variables = pr, tasmax, tasmin
        SDmethods = [
            SD : BCSD_4km*200[1-2].nc : 0 : 0,1
            CA : BCCA_4km*200[1-2].nc : 0 : 0,1
            SDe: BCSD_4km*200[1-2].nc : 0 : 0,1
            CAe: BCCA_4km*200[1-2].nc : 0 : 0,1
            ]
        gridsMatch=True
        lat = 34,44
        lon = -115,-100

    full_run.txt example:
        
    NOTE: format of SDmethods lines = 
        <DirectoryName>:<glob_search_pattern>:<dataset_indices>:<variable_indices>
        main patterns are : separated
        indices lists are , separated
        whitespace is stripped off all items
        Ends with a ] on its own line. 
        
AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    0.1

    
"""
import sys
import os
import traceback
import argparse
import re
import glob
import time

import numpy as np
import matplotlib.pyplot as plt

from bunch import Bunch
from nc_reader import NC_Reader
import swim_io as nc
import date_fun

debug=False

def default_short():
    '''Set up a default set of metadata to work on'''
    yrs="200[0-2]"
    obs=Bunch(dir="/Volumes/G-SAFE/usbr/statistical/obs/uw.0625/",search="*."+yrs+".nc")
    datasets=["narr","ncep"]
    varnames=["pr","tasmax","tasmin"]
    datasets=["narr","ncep"]
    SD=[Bunch(dir="SD",search="BCSD_4km*"+yrs+".nc",usedata=[0],usevars=[0,1]),
        Bunch(dir="CA",search="BCCA_4km*"+yrs+".nc",usedata=[0],usevars=[0,1]),
        Bunch(dir="SDe",search="BCSD_4km*"+yrs+".nc",usedata=[0],usevars=[0,1]),
        Bunch(dir="CAe",search="BCCA_4km*"+yrs+".nc",usedata=[0],usevars=[0,1])]
    gridsMatch=False
    lat=[34.5,43.5]
    lon=[-111.5,-100.5]
    return Bunch(obs=obs,datasets=datasets,varnames=varnames,
                SD=SD,gridsMatch=gridsMatch,geo=[lat[0],lat[1],lon[0],lon[1]])

def default_input():
    '''Set up a default set of metadata to work on'''
    return default_short()
    yrs="200[0-8]"
    obs=Bunch(dir="/Volumes/G-SAFE/usbr/statistical/obs/uw.0625/",search="*."+yrs+".nc")
    datasets=["narr","ncep"]
    varnames=["pr","tasmax","tasmin"]
    datasets=["narr","ncep"]
    SD=[Bunch(dir="SD",search="BCSD_4km*"+yrs+".nc",usedata=[0,1],usevars=[0,1,2]),
        Bunch(dir="CA",search="BCCA_4km*"+yrs+".nc",usedata=[0,1],usevars=[0,1,2]),
        Bunch(dir="SDe",search="BCSD_4km*"+yrs+".nc",usedata=[0,1],usevars=[0,1,2]),
        Bunch(dir="CAe",search="BCCA_4km*"+yrs+".nc",usedata=[0,1],usevars=[0,1,2]),
        Bunch(dir="SDe0",search="BCSD_4km*"+yrs+".nc",usedata=[0,1],usevars=[0,1,2]),
        Bunch(dir="CAe0",search="BCCA_4km*"+yrs+".nc",usedata=[0,1],usevars=[0,1,2]),
        Bunch(dir="SDe1",search="BCSD_4km*"+yrs+".nc",usedata=[0,1],usevars=[0,1,2]),
        Bunch(dir="CAe1",search="BCCA_4km*"+yrs+".nc",usedata=[0,1],usevars=[0,1,2])]
    gridsMatch=False
    lat=[34.5,43.5]
    lon=[-111.5,-100.5]
    return Bunch(obs=obs,datasets=datasets,varnames=varnames,
                SD=SD,gridsMatch=gridsMatch,geo=[lat[0],lat[1],lon[0],lon[1]])

def parse_inputfile(filename):
    '''Read which data to work on from an input file'''
    f=open(filename,'ru')
    defaults=default_input()
    SD=defaults.SD
    obs=defaults.obs
    datasets=defaults.datasets
    varnames=defaults.varnames
    gridsMatch=defaults.gridsMatch
    geo=defaults.geo
    for l in f:
        if re.match("\s*obs.*",l):
            data=l.split('=')[1].split(':')
            obs=Bunch(dir=data[0].strip(),search=data[1].strip())
        elif re.match("\s*datasets.*",l):
            datasets=[d.strip() for d in l.split('=')[1].split(',')]
        elif re.match("\s*variables.*",l):
            varnames=[d.strip() for d in l.split('=')[1].split(',')]
        elif re.match("\s*SDmethods\s*=.*",l):
            SD=list()
            l=f.next()
            while not re.match("\s*]\s*",l):
                data=l.split(':')
                l=f.next()
                SD.append(Bunch(dir=data[0].strip(),
                                search=data[1].strip(),
                                usedata=np.cast['i'](data[2].strip().split(',')),
                                usevars=np.cast['i'](data[3].strip().split(',')) ))
        elif re.match("\s*gridsMatch\s*=.*",l):
            gridsMatch=bool(l.split('=')[1].strip())
        elif re.match("\s*lat\s*=.*",l):
            lat=np.cast['f'](l.split('=')[1].split(','))
        elif re.match("\s*lon\s*=.*",l):
            lon=np.cast['f'](l.split('=')[1].split(','))
        
    f.close()
    return Bunch(obs=obs,datasets=datasets,varnames=varnames,
                SD=SD,gridsMatch=gridsMatch,geo=[lat[0],lat[1],lon[0],lon[1]])

def setup_datareaders(info):
    '''Sets up a structure of nc_reader objects to iterate over.'''
    # first set up the obs/verification data readers on for each possible variable
    obs_data=list()
    for var in info.varnames:
        print(info.obs.dir+'/'+var+'/'+info.obs.search,var)
        obs_data.append(NC_Reader(info.obs.dir+'/'+var+'/'+info.obs.search,readvars=[var],
                        geo_subset=info.geo,ntimes=365,latvar="lat",lonvar="lon"))
    size=obs_data[0].x.shape

    # then set up a list of lists of lists for the Statistical Downscaling readers
    # using lists so that each one can be a different length so we can do different combinations
    #  of variables and datasets for each method
    # top level are the different stat methods (CA,SD,SDe,...)
    sd_data=list()
    for current in info.SD:
        curdata=list()
        # second level are the different datasets (ncep,narr)
        for dataset in current.usedata:
            vardata=list()
            # third level are the different variables (pr,tasmin,tasmax)
            for usethis in current.usevars:
                var=info.varnames[usethis]
                # because this step takes a while, print out the current 
                # variable we are working on to make the user feel better
                print(current.dir,info.datasets[dataset],var)
                # add the current reader to the variable level list
                if info.gridsMatch:
                    # print(info.geo)
                    # checkfiles=glob.glob(current.dir+'/'+info.datasets[dataset]+'/'+var+'/'+current.search)
                    # print(current.dir+'/'+info.datasets[dataset]+'/'+var+'/'+current.search)
                    # lat=nc.read_nc(checkfiles[0],"lat").data
                    # lon=nc.read_nc(checkfiles[0],"lon").data
                    # print(lat.min(),lat.max(),lon.min()-360,lon.max()-360)
                    vardata.append(NC_Reader(current.dir+'/'+info.datasets[dataset]+'/'+var+'/'+current.search,readvars=[var],
                                    geo_subset=info.geo,ntimes=365,latvar="lat",lonvar="lon"))
                else:
                # the current reader will match this dataset with the geographic information in the
                # obs dataset if the grids don't match to begin with (info.gridsMatch).
                    vardata.append(NC_Reader(current.dir+'/'+info.datasets[dataset]+'/'+var+'/'+current.search,readvars=[var],
                                    geomatch_file=obs_data[0]._filenames[0],glatvar="lat",glonvar="lon",
                                    geo_subset=info.geo,ntimes=365,latvar="lat",lonvar="lon"))
            # add the current variable level list to the dataset level list
            curdata.append(vardata)
        # finally add the current dataset list to the SD method list
        sd_data.append(curdata)
                
    # return the obs list and the SD list of list of lists along with the size of the data array
    return(obs_data,sd_data,size)

def calc_errors(obslist,statlist,info,size):
    '''Calculate errors between all stat methods and observations'''
    
    # setup useful variables
    nSD=len(statlist)
    nvars=len(info.varnames)
    ndatasets=len(info.datasets)
    nyears=len(obslist[0]._filenames)
    ntimesteps=366*nyears #(that is actually slightly longer than necessary)
    points=np.array([[83,196],
                     [80,190],
                     [135,294],
                     [143,294],
                     [26,336],
                     [132,219],
                     [194,159]])
    points=np.array([])
    npoints=points.shape[0]
    nbins=100
    # output arrays are nSD_methods x 2(bias,RMS) x ndatasets x nvariables x (map=nx x ny)|(series=ntimesteps)
    # print("   Initializing output variables:")
    # print("    "+"x".join([str(i) for i in [nSD,2,ndatasets,nvars,size[0],size[1]]])+
    #         " = "+str(nSD*2*ndatasets*nvars*size[0]*size[1]/1024.0/1024)+" MB")
    # print("    "+"x".join([str(i) for i in [nSD,2,ndatasets,nvars,ntimesteps]])+
    #         " = "+str(nSD*2*ndatasets*nvars*ntimesteps/1024.0/1024)+" MB")
    map_out=np.zeros((nSD,2,ndatasets,nvars,size[0],size[1]),dtype=np.float64)
    time_out=np.zeros((nSD,2,ndatasets,nvars,ntimesteps))
    hist_out=np.zeros((nSD,ndatasets,nvars,npoints,nbins))
    obs_hist_out=np.zeros((nvars,npoints,nbins))
    
    hist_tmin=-35.
    hist_tmax=50
    hist_pmax=np.log10(50.)
    # these are just indices into the output arrays
    bias=0
    rms=1
    # counters for output (both screen and data)
    timestep=0
    outputinc=ntimesteps/20.0
    curoutput=outputinc
    # step through the observed dataset
    for pr in obslist[0]:
        tmax=obslist[1].next()
        tmin=obslist[2].next()
        # histogram calculations for obs (could be removed?)
        if (obslist[0]._curfile<len(obslist[0]._filenames)
            and re.match(".*2003.nc",obslist[0]._filenames[obslist[0]._curfile])  # skip 2003
            and (obslist[0].posinfile>210) and obslist[0].posinfile<244):         # in august
            pass
        else:
            for curpoint in range(npoints):
                if pr[0][points[curpoint,0],points[curpoint,0]] >= 1:
                    curval=np.round(np.log10(pr[0][points[curpoint,0],points[curpoint,1]])/hist_pmax*nbins)
                    if curval<0:curval=0
                    if curval>=nbins:
                        curval=nbins-1
                    obs_hist_out[0,curpoint,curval]+=1
            
                curval=np.round((tmax[0][points[curpoint,0],points[curpoint,1]]-hist_tmin)/(hist_tmax-hist_tmin)*nbins)
                if curval<0:curval=0
                if curval>=nbins:curval=nbins-1
                obs_hist_out[1,curpoint,curval]+=1

                curval=np.round((tmin[0][points[curpoint,0],points[curpoint,1]]-hist_tmin)/(hist_tmax-hist_tmin)*nbins)
                if curval<0:curval=0
                if curval>=nbins:curval=nbins-1
                obs_hist_out[2,curpoint,curval]+=1
            
        # status output to screen
        if timestep>=curoutput:
            print(np.round(float(timestep)/ntimesteps*100),end="% ")
            sys.stdout.flush()
            curoutput+=outputinc
        # create a list for the observed data
        obs=[pr[0],tmax[0],tmin[0]]
        # now step through the SD methods
        for i in range(nSD):
            statd=0
            # step through the datasets for the current method
            for d in info.SD[i].usedata:
                statvar=0
                # finally step through the variables for the current method.dataset
                for v in info.SD[i].usevars:
                    try:
                        statdata=statlist[i][statd][statvar].next()[0]
                        statvar+=1 #incase the SD data set doesn't use every var
                        if (obslist[0]._curfile<len(obslist[0]._filenames)
                            and re.match(".*2003.nc",obslist[0]._filenames[obslist[0]._curfile]) # skip 2003 (problems in)
                            and (obslist[0].posinfile>210) and obslist[0].posinfile<244):        # in august (uw dataset)
                            pass
                            # print(obslist[0]._filenames[obslist[0]._curfile])
                        else:
                            # print(statdata.shape,obs[v].shape)
                            errors=statdata-obs[v]
                            bad_data=np.abs(errors)>1000
                            errors[bad_data]=0
                            errors=np.ma.array(errors,mask=bad_data)
                            if info.varnames[v]=="pr":
                                gain=365.0
                            else:
                                gain=1.0
                            map_out[i,bias,d,v,...]+=errors *gain
                            map_out[i,rms,d,v,...]+=(errors*errors)# *gain
                            time_out[i,bias,d,v,timestep]=errors.mean()
                            time_out[i,rms,d,v,timestep]=np.sqrt((errors*errors).mean())
                            
                            # histogram calculations (could be removed?)
                            for curpoint in range(npoints):
                                if info.varnames[v]=="pr":
                                    if statdata[points[curpoint,0],points[curpoint,1]] >= 1:
                                        curval=np.round(np.log10(statdata[points[curpoint,0],points[curpoint,1]])/hist_pmax*nbins)
                                        if curval<0:curval=0
                                        if curval>=nbins:
                                            curval=nbins-1
                                        hist_out[i,d,v,curpoint,curval]+=1
                                else:
                                    curval=np.round((statdata[points[curpoint,0],points[curpoint,1]]-hist_tmin)/(hist_tmax-hist_tmin)*nbins)
                                    if curval<0:curval=0
                                    if curval>=nbins:curval=nbins-1
                                    hist_out[i,d,v,curpoint,curval]+=1
                                    
                            if debug and info.SD[i].dir=="SD":
                                plt.clf()
                                plt.subplot(221)
                                plt.imshow(statdata,origin="lower")
                                plt.title("Statistical Downscaling:  "+info.SD[i].dir)
                                plt.colorbar()
                                plt.subplot(222)
                                plt.imshow(obs[v],origin="lower")
                                plt.title("Observations")
                                plt.colorbar()
                                plt.subplot(223)
                                plt.imshow(errors,origin="lower")
                                plt.title("RMS:"+str(time_out[i,rms,d,v,timestep])[0:5]+"  Bias:"+str(time_out[i,bias,d,v,timestep])[0:5])
                                plt.colorbar()
                                plt.draw()
                                plt.show()
                                junk=raw_input()
                    except StopIteration:
                        print("Ran out of stat data:")
                        print(v,i,statd,timestep)
                statd+=1 #incase the SD data set doesn't use every dataset
        if (obslist[0]._curfile<len(obslist[0]._filenames)
            and re.match(".*2003.nc",obslist[0]._filenames[obslist[0]._curfile]) 
            and (obslist[0].posinfile>210) and obslist[0].posinfile<244):
            pass
            # print(obslist[0]._filenames[obslist[0]._curfile])
        else:
            timestep+=1
    time_out=time_out[:,:,:,:,:timestep]
    map_out/=timestep
    map_out[:,rms,...]=np.sqrt(map_out[:,rms,...])
    print("")
    return (map_out,time_out,obs_hist_out,hist_out)

def write_map(mapdata,filename,datatype="rms"):
    '''Write output map to netcdf file'''
    try:
        nc.write(filename,mapdata)
    except:
        print("IO error writing:"+str(filename+".nc"))
        
    if datatype=="rms":
        vmin=0
        tmp=mapdata[mapdata<1000].copy()
        tmp.sort()
        vmax=tmp[len(tmp)*0.99]
    else:
        tmp=mapdata[mapdata<1000].copy()
        tmp.sort()
        vmax=tmp[len(tmp)*0.99]
        vmin=tmp[len(tmp)*0.01]
        if vmax>np.abs(vmin):
            vmin=-vmax
        else:
            vmax=-vmin
    mask=np.abs(mapdata)<vmax*0.05
    mapdata=np.ma.array(mapdata,mask=mask)
    plt.clf()
    plt.imshow(mapdata,origin="lower",vmin=vmin,vmax=vmax)
    plt.colorbar()
    plt.title(filename)
    plt.savefig(filename+'.png')

def findfirstleap(info):
    files=glob.glob(info.obs.dir+info.obs.search)
    files.sort()
    for i,f in enumerate(files):
        year=int(f.split('.')[-2])
        if (year%4)==0:
            return i
    return -1

def conv2monthly(data,info):
    n=len(data)
    nyears=np.round(n/365.25)
    monthsize=np.array([31,28,31,30,31,30,31,31,30,31,30,31]).repeat(nyears+1)
    output=np.zeros(len(monthsize))
    firstleap=findfirstleap(info)
    if firstleap != -1:
        leapmonths=(np.arange(nyears/4.)+firstleap)*12+1
        leapmonths=leapmonths[where(leapmonths<len(monthsize))]
        monthsize[leapmonths]+=1
    curmonth=0
    nextmonth=monthsize[curmonth]
    for i in range(n):
        if i>nextmonth:
            output[curmonth]/=monthsize[curmonth]
            curmonth+=1
            nextmonth+=monthsize[curmonth]
        output[curmonth]+=data[i]
    return(output[:curmonth])
    

def write_series(seriesdata,filename,info,datatype="rms",starttime=None):
    '''Write output time series to netcdf file'''
    try:
        nc.write(filename,seriesdata)
    except:
        print("IO error writing:"+str(filename+".nc"))
    
    plt.clf()
    ndays=seriesdata.shape[0]
    if datatype=="rms":
        vmin=0
        tmp=seriesdata[seriesdata<1000].copy()
        tmp.sort()
        vmax=tmp[len(tmp)*0.98]
    else:
        tmp=seriesdata[seriesdata<1000].copy()
        tmp.sort()
        vmax=tmp[len(tmp)*0.98]
        vmin=tmp[len(tmp)*0.02]
        if vmax>np.abs(vmin):
            vmin=-vmax
        else:
            vmax=-vmin
    if starttime:
        xvals=date_fun.mjd2datetime(np.arange(ndays)+starttime)
    else:
        xvals=np.arange(ndays)
    monthly=conv2monthly(seriesdata,info)
    plt.subplot('211')
    plt.plot(xvals,seriesdata,'.')
    plt.plot([xvals[0],xvals[-1]],[0,0])
    fig=plt.gcf()
    fig.autofmt_xdate()
    # plt.xlim(0,len(seriesdata))
    plt.ylim(vmin,vmax)
    plt.title(filename+" Error="+str(seriesdata.mean())[0:5])
    plt.subplot('212')
    xvals=(np.arange(len(monthly))+0.5)/12.0
    plt.plot(xvals,monthly)
    plt.plot([xvals[0],xvals[-1]],[0,0])
    # plt.xlim(xvals[0],xvals[-1])
    plt.savefig(filename+'.png')

def write_hist(obs_hist,sd_hist,filename,vmin=0.0,vmax=None,xlog=False):
    nbins=obs_hist.shape[1]
    if vmax==None:
        vmax=nbins
    stepsize=(vmax-vmin)/nbins
    xvals=np.arange(nbins)/nbins*(vmax-vmin)+vmin
    xvals=np.arange(vmin,vmax,stepsize)
    if xlog:
        xvals=10**xvals
    points=np.array([[83,196],
                     [80,190],
                     [135,294],
                     [143,294],
                     [26,336],
                     [132,219],
                     [194,159]])
    try:
        nc.write(filename+"_hists",sd_hist)
    except:
        print("IO error writing:"+str(filename+"_hists.nc"))
    npoints=obs_hist.shape[0]
    for point in range(npoints):
        plt.clf()
        ax=plt.subplot(111)
        ax.plot(xvals,obs_hist[point,:],color="b",linewidth=2,label="Observations")
        ax.plot(xvals,sd_hist[point,:],color="r",linewidth=2,label="Downscaled")
        plt.xlim(xvals[0],xvals[-1])
        if xlog:
            ax.set_xscale('log')
        if obs_hist[point,0]>obs_hist[point,1:].max():
            newmax=max(obs_hist[point,1:].max(),sd_hist[point,1:].max())
            plt.ylim(0,newmax*1.1)
        plt.legend()
        x=str(points[point,1])
        y=str(points[point,0])
        plt.title(filename.replace("_"," ")+"  x="+x+"  y="+y)
        plt.draw()
        plt.savefig(filename+"_"+x+"_"+y+".png")

def write_results(results,info,starttime=None):
    '''Iterate over results and output map and time series results'''
    maps=results[0]
    series=results[1]
    obs_hist=results[2]
    # nc.write("obs_hists",obs_hist)
    sd_hist=results[3]
    bias=0
    rms=1
    hist_tmin=-35.
    hist_tmax=50
    hist_pmax=np.log10(50.)
    nSD=maps.shape[0]
    for i in range(nSD):
        for d in info.SD[i].usedata:
            dataset=info.datasets[d]
            for v in info.SD[i].usevars:
                curvar=info.varnames[v]
                write_map(maps[i,bias,d,v,...],'_'.join([info.SD[i].dir,dataset,curvar,"bias","map"]),datatype="bias")
                write_map(maps[i,rms,d,v,...],'_'.join([info.SD[i].dir,dataset,curvar,"RMS","map"]),datatype="rms")
                write_series(series[i,bias,d,v,...],'_'.join([info.SD[i].dir,dataset,curvar,"bias","series"]),info,datatype="bias",starttime=starttime)
                write_series(series[i,rms,d,v,...],'_'.join([info.SD[i].dir,dataset,curvar,"RMS","series"]),info,datatype="rms",starttime=starttime)
                # if curvar=="pr":
                #     write_hist(obs_hist[v,...],sd_hist[i,d,v,...],'_'.join([info.SD[i].dir,dataset,curvar]),vmin=0,vmax=hist_pmax,xlog=True)
                # else:
                #     write_hist(obs_hist[v,...],sd_hist[i,d,v,...],'_'.join([info.SD[i].dir,dataset,curvar]),vmin=hist_tmin,vmax=hist_tmax)
                    

def main (filename=None):

    starttime=date_fun.date2mjd(2000,1,1,0,0)
    if filename:
        info=parse_inputfile(filename)
    else:
        info=default_input()
    print("Run Info:")
    print(info)
    print("Setting up input file readers")
    obs,stat,size=setup_datareaders(info)
    print("Calculating statistics")
    results=calc_errors(obs,stat,info,size)
    print("Writing output results")
    write_results(results,info,starttime=starttime)

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='A simple statistical downscaling comparison. ')
        parser.add_argument('filename',action='store')
        parser.add_argument('-v', '--version',action='version',
                version='0.1')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        exit_code = main(args.filename)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print( 'ERROR, UNEXPECTED EXCEPTION')
        print( str(e))
        traceback.print_exc()
        os._exit(1)
