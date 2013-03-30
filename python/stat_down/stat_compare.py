#!/usr/bin/env python
from __future__ import print_function
import glob
import sys
import traceback
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from nc_reader import NC_Reader
from swim_io import NC_writer
import date_fun
import swim_io as mygis

wrf_geo='/glade/scratch/gutmann/usbr/hybrid/baseline_info/4km_wrf_output.nc'
# read in XLAT and XLONG
# in main wrf files, data var = "RAINNC" units are mm
# sd_geo='/Volumes/G-SAFE/usbr/statistical/down/CA/ncep/pr/BCCA_4km_T62_prate.sfc.gauss.2000.nc'
# read in lat and lon, lon=lon-360, (lat,lon)=np.meshgrid(lat,lon)
# wrf_subset = sd_geo.data[216:451,242:622]
# in main SD/CA files data var = "pr" units=(mm?)
# wrffiles=glob.glob("/Volumes/G-SAFE/usbr/wrf4km_daily_precip/NARR*.nc")
def output_movie(file_prefix,data,dataareprecip=False):
    if dataareprecip:
        # test=data[0:365,:,:]
        # good=(test<100)&(test>0)
        # mn=test[good].min()
        # mx=test[good].max()
        mn=0
        mx=50
    else:
        winter=data[0:40,:,:]
        summer=data[180:220,:,:]
    
        wintermin=winter[winter<100].min()
        wintermax=winter[winter<100].max()
        summermin=summer[summer<100].min()
        summermax=summer[summer<100].max()
        minrange=(summermin-wintermin)/2
        maxrange=(summermax-wintermax)/2
    
    outputdir=file_prefix+"_movie/"
    try:
        os.mkdir(outputdir)
    except:
        pass
    lat=np.array([30,45])
    lon=np.array([-120,-95])
    
    for i in range(data.shape[0]):
        print(str(i)+"/"+str(data.shape[0]))
        if i<10:
            iprefix="000"
        elif i<100:
            iprefix="00"
        elif i<1000:
            iprefix="0"
        else:
            iprefix=""
        curfile=file_prefix.replace('_',' ')+' '+iprefix+str(i)
        if not dataareprecip:
            curpoint= 1-np.cos((i/365.25)*2*np.pi)
            mn=curpoint*minrange+wintermin
            mx=curpoint*maxrange+wintermax
        map_vis(data[i,:,:],lat,lon,title=curfile,vmin=mn,vmax=mx,outputdir=outputdir)

def runfor(SD_data,precdiff=False,movie=None,dataareprecip=False,dataaretmin=False, usetimes=None):
    # run through the first iteration to set up all the variables
    if precdiff:
        if SD_data.posinfile>0:
            SD_data.posinfile-=1
            lastdata=SD_data.next()[0]
            curdata=SD_data.next()[0]-lastdata
        else:
            curdata=SD_data.next()[0]
    else:
        curdata=SD_data.next()[0]
    n=1.0
    ny,nx=curdata.shape
    print(curdata.shape)
    ntimes=len(SD_data._filenames)*366
    if len(SD_data._filenames)>50:
        ntimes/=12
    
    wet_threshold=2.54
    # wet_threshold=0.1
    wet_threshold=0.0
    ndays_sum=5
    
    alldata=np.zeros((ntimes,ny,nx),dtype=np.float32)
    alldata[n-1,:,:]=curdata
    sddata=curdata.copy()
    wetdata=np.choose(sddata>=wet_threshold,(0,1))
    totalwet=np.zeros((ny,nx),dtype=np.float32)
    totaldry=np.zeros((ny,nx),dtype=np.float32)
    lastwet=np.zeros((ny,nx),dtype=np.float32)
    curdry_spell=np.zeros((ny,nx),dtype=np.float32)
    maxdry_spell=np.zeros((ny,nx),dtype=np.float32)
    nwet=np.zeros(wetdata.shape)
    nwet[np.where(wetdata==1)]=1
    ndry=np.zeros(wetdata.shape)
    ndry[np.where(wetdata==0)]=1
    nwet2=np.zeros(wetdata.shape)
    nwet2[np.where(wetdata==1)]=1
    lastdata=curdata

    curtime=0
    # calc stats from every timestep in the file
    for curdata in SD_data:
        curtime+=1
        if (usetimes==None) or (curtime in usetimes):
            curdata=curdata[0]
            # check for missing data, this may not be necessary anymore, but originally there were some problems. 
            if curdata[10,10]<3000:
                if (n%100)==0:
                    print(str(n)+" ",end="")
                    sys.stdout.flush()
                n+=1
                if n>=ntimes:
                    ntimes*=2
                    new_alldata=np.zeros((ntimes,ny,nx),dtype=np.float32)
                    new_alldata[:n,...]=alldata[:n,...]
                    del(alldata)
                    alldata=new_alldata
                if precdiff:
                    alldata[n-1,:,:]=curdata-lastdata
                    sddata+=curdata-lastdata
                    wetdata=np.choose((curdata-lastdata)>wet_threshold,(0,1))
                    lastdata=curdata
                else:
                    alldata[n-1,:,:]=curdata
                    sddata+=curdata
                    if dataareprecip:
                        wetdata=np.choose(curdata>wet_threshold,(0,1))
                    else:
                        wetdata=np.choose(curdata>=0,(0,1))

                # Average wet spell length calculation:
                # add one to the total number of wet days everywhere it rained
                totalwet+=wetdata
                tmp=np.where((lastwet==0) & (wetdata==1))
                # if this is a new wet period (yesterday was dry) add one to the wet period count
                nwet[tmp]+=1
                tmpwet=np.where(wetdata==1)
                nwet2[tmpwet]+=1

                # Average dry spell length calculation: (inverse of the above calculation)
                # add one to the total number of dry days everywhere it was dry
                totaldry+=1-wetdata
                tmp=np.where((lastwet==1) & (wetdata==0))
                # if this is a new dry period (yesterday was wet) add one to the dry period count
                ndry[tmp]+=1
                
                lastwet=wetdata
                
                # maximum dry spell length calculation:
                # first assume everywhere is dry and add one to the dry spell length
                curdry_spell+=1
                # then find places where it did rain and set dry spell length to 0
                tmp=np.where(wetdata>0)
                if len(tmp[0])>0:
                    curdry_spell[tmp]=0
                # finally if anywhere has a longer dryspell than the longest on record
                #  make the current value the new record for those locations
                tmp=np.where(curdry_spell>maxdry_spell)
                if len(tmp[0])>0:
                    maxdry_spell[tmp]=curdry_spell[tmp]
        else:
            # the current time slice is not in the correct month(s)
            pass
            # technically we should be able to increment SDdata.posinfile
            # but that is a little tricky given >1 file 
            # (but might be worth it if this takes a while)
    
    print("ndry_periods",ndry.mean())     
    print("ndry_length",totaldry.mean())     
    print(maxdry_spell.max())
    print(maxdry_spell.min())
    print("Got All Data")
    alldata=alldata[0:n,:,:]
    sz=alldata.shape
    ndaydata=np.zeros((sz[0],sz[1],sz[2]))
    halfsum=np.floor(ndays_sum/2)
    isodd=ndays_sum%2
    for i in range(ndays_sum):
        print(ndaydata[halfsum:-(halfsum+isodd),:,:].shape, alldata[i:i-ndays_sum,:,:].shape)
        if i==ndays_sum-1:
            ndaydata[halfsum:-(halfsum+isodd),:,:]+=alldata[i+1:,:,:]
        else:
            ndaydata[halfsum:-(halfsum+isodd),:,:]+=alldata[i:i-ndays_sum,:,:]
    if movie:
        print("making the movie...")
        output_movie(movie,alldata[-365*4-1:,:,:],dataareprecip=dataareprecip)
    print("Sorting")
    alldata.sort(axis=0)
    # ndaydata.sort(axis=0)
    print("...sorted")
    # print(ndaydata.max())
    # calculate min or max 10th/90th, 5th/95th, and 1st/99th percentiles
    if dataaretmin:
        maxprec90=alldata[np.int(n*0.10),:,:]
        maxprec95=alldata[np.int(n*0.05),:,:]
        maxprec99=alldata[np.int(n*0.01),:,:]
        for test in range(np.int((n-ndays_sum)*0.01)):
            curmaxpos=np.argmax(ndaydata,axis=0)
            for eliminate in range(ndays_sum):
                ndaydata[curmaxpos+int(eliminate-halfsum),:,:]=0
        maxprec99_nday=np.max(ndaydata,axis=0)
    else:
        maxprec90=alldata[np.int(n*0.90),:,:]
        maxprec95=alldata[np.int(n*0.95),:,:]
        maxprec99=alldata[np.int(n*0.99),:,:]
        for test in range(np.int((n/ndays_sum)*0.01)):
            curmaxpos=np.argmax(ndaydata,axis=0)
            for eliminate in range(ndays_sum):
                ndaydata[curmaxpos+int(eliminate-halfsum),:,:]=0
            print(test,ndaydata.max())
            
        maxprec99_nday=np.max(ndaydata,axis=0)
        
    print("calculating stddev")
    stddev=maxprec90.copy()
    # stddev=np.std(alldata,axis=0,ddof=1)
    # wrap up the stats (mostly normalizing by n)
    print("calculating wet day statistics")
    nwet[np.where(nwet==0)]=1
    ndry[np.where(ndry==0)]=1
    wet_duration=totalwet/nwet.astype('f')
    dry_duration=totaldry/ndry.astype('f')
    print("done with stats")
    return(sddata/n,nwet2/n,wet_duration,stddev,maxprec90,maxprec95,maxprec99,maxprec99_nday,maxdry_spell,dry_duration)


def map_vis(data,lat,lon,title="",vmin=None,vmax=None,dx=4000.0,barlabel=None,outputdir=None,geo=None):
    print("visualizing : "+title)
    if (vmin < -30):
        vmin=-30
    if len(data.shape)>2:
        plotdata=data[0,:,:]
    else:
        plotdata=data
    plt.clf()
    ny,nx=plotdata.shape
    if ny>nx:
        plotdata=plotdata.T
        ny,nx=plotdata.shape
    if re.match(".*async.*",title):
        m = Basemap(width=nx*dx,height=ny*dx,
                    rsphere=(6378137.00,6356752.3142),\
                    resolution='l',area_thresh=1000.,projection='lcc',\
                    lat_1=33.,lat_2=45.,lat_0=39.,lon_0=-105.875)
    else:
        m = Basemap(width=nx*dx,height=ny*dx,
                    rsphere=(6378137.00,6356752.3142),\
                    resolution='l',area_thresh=1000.,projection='lcc',\
                    lat_1=33.,lat_2=45.,lat_0=39,lon_0=-107.)
    if geo!=None:
        m = Basemap(projection='merc',llcrnrlat=geo[0],urcrnrlat=geo[1],\
                    llcrnrlon=geo[2],urcrnrlon=geo[3],lat_ts=(geo[0]+geo[1])/2)
        # 
        # m = Basemap(width=nx*dx,height=ny*dx,
        #             rsphere=(6378137.00,6356752.3142),\
        #             resolution='l',area_thresh=1000.,projection='lcc',\
        #             lat_1=33.,lat_2=45.,lat_0=39,lon_0=-107.)
        
    mapimg=m.imshow(plotdata,vmin=vmin,vmax=vmax)
    m.drawparallels(np.arange(lat.min(),lat.max(),2.),labels=[1,1,0,0],dashes=[1,4])
    m.drawmeridians(np.arange(lon.min(),lon.max(),4.),labels=[0,0,0,1],dashes=[1,4])
    m.drawstates(linewidth=1.5)
    cbar=m.colorbar()
    if barlabel:
        cbar.set_label(barlabel)
    plt.title(title)
    if outputdir:
        plt.savefig(outputdir+title.replace(' ','_')+'.png')
    else:
        plt.savefig(title.replace(' ','_')+'.png')


def write_data(outprefix,sddata,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration,reverse_index=False,dataareprecip=False,dataaretmin=False,geo=None):
    if reverse_index:
        sddata=sddata.T
        wet_frac=wet_frac.T
        wet_duration=wet_duration.T
        maxdry=maxdry.T
    if m90[10,10]>100:
        sddata-=273.15
        m90-=273.15
        m95-=273.15
        m99-=273.15
    (Nx,Ny)=sddata.shape
    lat=np.array([30,45])
    lon=np.array([-120,-95])
    
    # write output data
    ncout=NC_writer(outprefix+'_mean_precip',Ny,Nx,var='precip')
    ncout.appendToVar(sddata)
    ncout.close()
    if dataareprecip:
        map_vis(sddata*365,lat,lon,outprefix.replace('_',' ')+" MAP",vmin=0,vmax=5*365.0,geo=geo)
    else:
        if dataaretmin:
            prefix="Min"
            mn=-13
            mx=10
        else:
            prefix="Max"
            mx=28
            mn=2
        # mx=sddata[sddata<100].max()
        # mn=sddata[sddata>-100].min()
        map_vis(sddata,lat,lon,outprefix.replace('_',' ')+" Mean "+prefix+" Temp.",vmin=mn,vmax=mx,geo=geo)
        
    ncout=NC_writer(outprefix+'_wet_frac',Ny,Nx,var='wetfrac')
    ncout.appendToVar(wet_frac)
    ncout.close()
    if dataareprecip:
        map_vis(wet_frac,lat,lon,outprefix.replace('_',' ')+" Wet Fraction",vmin=0.2,vmax=0.9,geo=geo)
        # map_vis(wet_frac,lat,lon,outprefix.replace('_',' ')+" Wet Fraction",vmin=0.6,vmax=1,geo=geo)
    else:
        try:
            mx=wet_frac[wet_frac<100].max()
            mn=wet_frac[wet_frac>-100].min()
        except:
            mx=wet_frac.max()
            mn=wet_frac.min()
            print("ERROR in max/min calc for wet_frac: ")
            print("  max="+str(mx)+"   min="+str(mn))
        map_vis(wet_frac,lat,lon,outprefix.replace('_',' ')+" Fraction "+prefix+" T above Freezing",vmin=mn,vmax=mx,geo=geo)
    
    ncout=NC_writer(outprefix+'_wet_duration',Ny,Nx,var='duration')
    ncout.appendToVar(wet_duration)
    ncout.close()
    if dataareprecip:
        map_vis(wet_duration,lat,lon,outprefix.replace('_',' ')+" Wet Day Duration",vmin=1,vmax=10,geo=geo)
    # else:
    #     map_vis(wet_duration,lat,lon,outprefix.replace('_',' ')+" Days in a row >0",vmin=1,vmax=4.5)
    
    ncout=NC_writer(outprefix+'_std_dev',Ny,Nx,var='stddev')
    ncout.appendToVar(stddev)
    ncout.close()
    if dataareprecip:
        map_vis(stddev,lat,lon,outprefix.replace('_',' ')+" Standard Deviation",vmin=0,vmax=10,geo=geo)
    else:
        # mx=stddev[stddev<100].max()
        # mn=stddev[stddev<100].min()
        mx=13.3
        mn=7.8
            
        map_vis(stddev,lat,lon,outprefix.replace('_',' ')+" Standard Deviation",vmin=mn,vmax=mx,geo=geo)
                
    
    ncout=NC_writer(outprefix+'_90percentile',Ny,Nx,var='precip')
    ncout.appendToVar(m90)
    ncout.close()
    if dataareprecip:
        map_vis(m90,lat,lon,outprefix.replace('_',' ')+" 90th percentile",vmin=0,vmax=15,geo=geo)
    else:
        # mx=m90[m90<100].max()
        # mn=m90[m90<100].min()
        if dataaretmin:
            percentile=" 10th"
            mx=0
            mn=-30
        else:
            percentile=" 90th"
            mx=38
            mn=17
        map_vis(m90,lat,lon,outprefix.replace('_',' ')+percentile+" percentile",vmin=mn,vmax=mx,geo=geo)
    
    ncout=NC_writer(outprefix+'_95percentile',Ny,Nx,var='precip')
    ncout.appendToVar(m95)
    ncout.close()
    if dataareprecip:
        map_vis(m95,lat,lon,outprefix.replace('_',' ')+" 95th percentile",vmin=0,vmax=25,geo=geo)
    else:
        # mx=m95[m95<100].max()
        # mn=m95[m95<100].min()
        if dataaretmin:
            percentile="  5th"
            mx=-3
            mn=-30
        else:
            percentile=" 95th"
            mx=41
            mn=19
        map_vis(m95,lat,lon,outprefix.replace('_',' ')+percentile+" percentile",vmin=mn,vmax=mx,geo=geo)
        
    ncout=NC_writer(outprefix+'_99percentile',Ny,Nx,var='precip')
    ncout.appendToVar(m99)
    ncout.close()
    if dataareprecip:
        map_vis(m99,lat,lon,outprefix.replace('_',' ')+" 99th percentile",vmin=0,vmax=35,geo=geo)
    else:
        # mx=m99[m99<100].max()
        # mn=m99[m99<100].min()
        if dataaretmin:
            percentile="  1st"
            mx=-5
            mn=-35
        else:
            percentile=" 99th"
            mx=43
            mn=21
        map_vis(m99,lat,lon,outprefix.replace('_',' ')+percentile+" percentile",vmin=mn,vmax=mx,geo=geo)

    ncout=NC_writer(outprefix+'_99_ndaypercentile',Ny,Nx,var='precip')
    ncout.appendToVar(maxprec99_nday)
    ncout.close()
    if dataareprecip:
        map_vis(maxprec99_nday,lat,lon,outprefix.replace('_',' ')+" 99th n-day percentile",vmin=0,vmax=45,geo=geo)
    else:
        # mx=m99[m99<100].max()
        # mn=m99[m99<100].min()
        if dataaretmin:
            percentile="  1st"
            mx=-5
            mn=-35
        else:
            percentile=" 99th"
            mx=43
            mn=21
        map_vis(maxprec99_nday/5,lat,lon,outprefix.replace('_',' ')+percentile+" n-day percentile",vmin=mn,vmax=mx,geo=geo)


    if dataareprecip:
        ncout=NC_writer(outprefix+'_dryspell',Ny,Nx,var='precip')
        ncout.appendToVar(maxdry)
        ncout.close()
        # bufferzone=30
        # mdrymin=maxdry[bufferzone:-bufferzone,bufferzone:-bufferzone].min()
        # mdrymax=maxdry[bufferzone:-bufferzone,bufferzone:-bufferzone].max()
        # if mdrymin<0:
        #     mdrymin=0
        mdrymin=8
        mdrymax=80
        map_vis(maxdry,lat,lon,outprefix.replace('_',' ')+" max dry spell",vmin=mdrymin,vmax=mdrymax,geo=geo)
        # if maxdry.max()>10:
        #     map_vis(maxdry,lat,lon,outprefix.replace('_',' ')+" max dry spell",vmin=0,vmax=90)
        # else:
        #     map_vis(maxdry,lat,lon,outprefix.replace('_',' ')+" max dry spell",vmin=0,vmax=8)
    ncout=NC_writer(outprefix+'_dry_duration',Ny,Nx,var='duration')
    ncout.appendToVar(dry_duration)
    ncout.close()
    if dataareprecip:
        map_vis(dry_duration,lat,lon,outprefix.replace('_',' ')+" Dry Day Duration",vmin=1,vmax=10,geo=geo)
    # else:
    #     map_vis(dry_duration,lat,lon,outprefix.replace('_',' ')+" Days in a row >0",vmin=1,vmax=4.5)
    

def obs2date(time):
    base_date=date_fun.date2mjd([1800,1,1,0,0])
    date_conversion=1.0
    dates=date_fun.mjd2date(time/time_conversion+basetime)
    return dates[1,:]

def SD2date(time):
    base_date=date_fun.date2mjd([1940,1,1,0,0])
    date_conversion=1.0/24
    dates=date_fun.mjd2date(time/time_conversion+basetime)
    return dates[1,:]
    
def wrf2date(time):
    try:
        ntimes=len(time)
    except TypeError:
        ntimes=1
        time=[time]
    months=np.zeros(ntimes)
    for i,t in enumerate(time):
        months[i]=int(str(t)[4:6])
    return months

def wrf4k_2date(time):
    try:
        ntimes=len(time)
    except TypeError:
        ntimes=1
        time=[time]
    months=np.zeros(ntimes)
    for i,t in enumerate(time):
        months[i]=int(t[5:7])
    return months

def make_usetimes(month,files,var="time",dateconverter=SD2date):
    '''Allows stat_compare to be run for individual months
    
    create a time series of indices to be used for the current month'''
    try:
        nmonths=len(month)
    except TypeError:
        month=[month] #make month iterable so we can check "date in month"
    timestamps=np.zeros((len(files)*366))
    i=0
    # loop over files reading the time variable and makeing a list of corresponding months
    for f in files:
        curtimes=mygis.read_nc(f,var).data
        dates=dateconverter(curtimes)
        # if we are past the length of the current array, double its size. 
        if (i+len(curtimes))>len(timestamps):
            sz*=2
            newtimestamps=np.zeros((sz))
            newtimestamps[:sz/2]=timestamps
            timestamps=newtimestamps
        timestamps[i:i+len(curtimes)]=dates
        i+=len(curtimes)
    timestamps=timestamps[:i]
    
    outputindices=np.zeros((len(timestamps)))
    n=0
    # now collect the list of times that are for the correct month
    for i,t in enumerate(timestamps):
        if t in month:
            outputindicies[n]=i
            n+=1
    return outputindicies[:n]


def main():
    skipobs=False
    runsd=True
    # stattype=['CA','SD','CAe','SDe','SDdry','SDwet','CAe0','SDe0','CAe1','SDe1',"CAe2","CAe3"]
    stattype=["SARe0","SARe1",'SDmon','CAe','SDe','CAe0','SDe0','CAe1','SDe1','SD','CA']
    stattype=['CAe0','SDe0','CAe1','SDe1','SD','CA']
    # stattype=['CAe','SDe','CAe0','SDe0','CAe1','SDe1','SD','CA']
    # stattype=['SDmon','CAe','SDe','CAe0','SDe0','CAe1','SDe1','SD','CA',
    #           'SDdry','SDwet','SDhot','SDcold','CAdry','CAwet','CAhot','CAcold']
    # stattype=['CAe','SDe','SDdry','SDwet','CAe1','SDe1']
    # stattype=['CAe','SDe','CAe0','SDe0','CAe1','SDe1']
    # stattype=['CAe1','SDe1','SARe1']
    # stattype=['SDe0']

    basedata=['ncep','narr']
    # basedata=['ncep']
    # basedata=['narr']

    # restype=['12km']
    # restype=['6km']
    # restype=['4km']
    restype=['12km','6km']

    vartype=['pr','tasmin','tasmax']
    # vartype=['pr']
    # vartype=['tasmin','tasmax']
        
    # makemovie=True
    # BCstatus=""
    BCstatus="BC"
    makemovie=False
    # lat=mygis.read_nc(wrf_geo,var="XLAT").data
    # lon=mygis.read_nc(wrf_geo,var="XLONG").data
    month=None
    monthname="" #note this could be e.g. DJF or January or JAN or "01" or ...
    usetimes=None
    test_SAR=True
    if test_SAR:
        # stattype=["SARe0","SARe1"]#,"SARdry","SARwet","SARhot","SARcold"]
        stattype=["SARe0"]
        restype=['12km','6km']
        basedata=['narr','ncep']
    latmin=34;latmax=43.9;lonmin=-115;lonmax=-100
    # latmin=35;latmax=37;lonmin=-111;lonmax=-109
    if runsd:
        for base in basedata:
            for stat in stattype:
                for res in restype:
                    for var in vartype:
                        print(base,stat,res,var)
                        if base is "ncep":
                            sdfilesearch="./"+stat+"/"+base+"/"+var+"/"+BCstatus+stat[:2]+"*_"+res+"*gauss*200[1-8]*.nc"
                        if base is "narr":
                            sdfilesearch="./"+stat+"/"+base+"/"+var+"/"+BCstatus+stat[:2]+"*_"+res+"*200[1-8]*.nc"
                        print(sdfilesearch)
                        outprefix="_".join(["base_stats",base,stat,res,var])
                        base_outprefix=outprefix
                    
                        checkfiles=glob.glob(sdfilesearch)
                        if month!=None:
                            dateconverter=SD2date
                            usetimes=make_usetimes(month,checkfiles,"time",dateconverter)
                            outprefix=base_outprefix+"_{}".format(monthname)
                        dataareprecip= (var=="pr")
                        if makemovie:
                            moviefilename=outprefix
                        else:
                            moviefilename=None
                        if len(checkfiles)>1:
                            try:
                                # SD_data=NC_Reader(sdfilesearch,checkfiles[0],wrf_geo,readvars=[var],
                                #                     ntimes=365,latvar="lat",lonvar="lon",subset=[10,-10,10,-10])
                                SD_data=NC_Reader(sdfilesearch,checkfiles[0],readvars=[var],ntimes=365,
                                                    latvar="lat",lonvar="lon",geo_subset=[latmin,latmax,lonmin,lonmax])
                                (meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)=runfor(SD_data,movie=moviefilename,
                                            dataareprecip=(var=="pr"),dataaretmin=(var=="tasmin"),usetimes=usetimes)
                                write_data(outprefix,meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration,
                                            dataareprecip=(var=="pr"),dataaretmin=(var=="tasmin"),geo=[latmin,latmax,lonmin,lonmax])
                            except Exception, e:
                                print('ERROR, UNEXPECTED EXCEPTION')
                                print(str(e))
                                traceback.print_exc()
                                print("Something went wrong with: "+outprefix.replace('_',' '))
                            
                        else:
                            print("NO FILES FOR : "+outprefix)
    if skipobs:
        return
    
    # finally, run for obs data
    dateconverter=obs2date
    for r in restype:
        for v in vartype:
            if r=="12km":
                obsfilesearch="../obs/maurer.125/"+v+"/*.200[1-8]*"
            elif r=="4km":
                obsfilesearch="../obs/uw.4km/"+v+"/*.200[1-8]*"
            elif r=="6km":
                obsfilesearch="../obs/uw.0625/"+v+"/*.200[1-8]*"
                
            checkfiles=glob.glob(obsfilesearch)
            if len(checkfiles)>0:
                # obs_data=NC_Reader(obsfilesearch,checkfiles[0],wrf_geo,
                #                     readvars=[v],ntimes=365,latvar="lat",lonvar="lon",subset=[10,-10,10,-10])
                obs_data=NC_Reader(obsfilesearch,checkfiles[0],readvars=[v],ntimes=365,
                                    latvar="lat",lonvar="lon",geo_subset=[latmin,latmax,lonmin,lonmax])
                (meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)=runfor(obs_data,
                                    dataareprecip=(v=="pr"),dataaretmin=(v=="tasmin"))
                outprefix="base_stats_obs_"+r+"_"+v
                write_data(outprefix,meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration,
                            dataareprecip=(v=="pr"),dataaretmin=(v=="tasmin"),geo=[latmin,latmax,lonmin,lonmax])
            else:
                print("No obs files for : "+str(v)+" "+r)
    
    # run for WRF data
    # if "4km" in restype:
    #     if "tasmax" in vartype:
    #         wrffilesearch="/Volumes/Data2/usbr/wrf_daily_tasmax*"
    #         checkfiles=glob.glob(wrffilesearch)
    #         if len(checkfiles)>0:
    #             wrf_data=NC_Reader(wrffilesearch,readvars=["tasmax"],
    #                                 ntimes=365,firstfile_timeinit=92,subset=[10,-10,10,-10])
    #             (meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)=runfor(wrf_data)
    #             outprefix="wrf_dyn_4km_tasmax"
    #             write_data(outprefix,meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration,reverse_index=True)
    #         else:
    #             print("No files for : wrf tasmax 4km")
    # 
    #     if "tasmin" in vartype:
    #         wrffilesearch="/Volumes/Data2/usbr/wrf_daily_tasmin*"
    #         checkfiles=glob.glob(wrffilesearch)
    #         if len(checkfiles)>0:
    #             wrf_data=NC_Reader(wrffilesearch,readvars=["tasmin"],
    #                                 ntimes=365,firstfile_timeinit=92,subset=[10,-10,10,-10])
    #             (meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)=runfor(wrf_data,dataaretmin=True)
    #             outprefix="wrf_dyn_4km_tasmin"
    #             write_data(outprefix,meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration,reverse_index=True,dataaretmin=True)
    #         else:
    #             print("No files for : wrf tasmin 4km")
    # 
    #     if "pr" in vartype:
    #         wrffilesearch="/Volumes/G-SAFE/usbr/wrf4km_daily_precip/NARR*"
    #         checkfiles=glob.glob(wrffilesearch)
    #         if len(checkfiles)>0:
    #             wrf_data=NC_Reader(wrffilesearch,readvars=["RAINNC"],
    #                                 ntimes=365,firstfile_timeinit=92,subset=[10,-10,10,-10])
    #             (meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)=runfor(wrf_data,precdiff=True,dataareprecip=True)
    #             outprefix="wrf_dyn_4km_precip"
    #             write_data(outprefix,meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration,dataareprecip=True,reverse_index=True)
    #         else:
    #             print("No files for : wrf precip 4km")
    # # run for wrf data regridded to 12km
    # if "12km" in restype:
    #     if "tasmax" in vartype:
    #         wrffilesearch="/Volumes/Data2/usbr/added*tasmax*.nc"
    #         checkfiles=glob.glob(wrffilesearch)
    #         if len(checkfiles)>0:
    #             wrf_data=NC_Reader(wrffilesearch,checkfiles[0],wrf_geo,readvars=["tasmax"],
    #                                 ntimes=365,firstfile_timeinit=92,latvar="lat",lonvar="lon",subset=[10,-10,10,-10])
    #             (meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)=runfor(wrf_data)
    #             outprefix="wrf_dyn_12km_tasmax"
    #             write_data(outprefix,meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)
    #         else:
    #             print("No files for : wrf tasmax 12km")
    #     
    #     if "tasmin" in vartype:
    #         wrffilesearch="/Volumes/Data2/usbr/added*tasmin*.nc"
    #         checkfiles=glob.glob(wrffilesearch)
    #         if len(checkfiles)>0:
    #             wrf_data=NC_Reader(wrffilesearch,checkfiles[0],wrf_geo,readvars=["tasmin"],
    #                                 ntimes=365,firstfile_timeinit=92,latvar="lat",lonvar="lon",subset=[10,-10,10,-10])
    #             (meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)=runfor(wrf_data,dataaretmin=True)
    #             outprefix="wrf_dyn_12km_tasmin"
    #             write_data(outprefix,meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration,dataaretmin=True)
    #         else:
    #             print("No files for : wrf tasmin 12km")
    #     
    #     if "pr" in vartype:
    #         wrffilesearch="/Volumes/G-SAFE/usbr/wrf4km_daily_precip/added*.nc"
    #         checkfiles=glob.glob(wrffilesearch)
    #         if len(checkfiles)>0:
    #             wrf_data=NC_Reader(wrffilesearch,checkfiles[0],wrf_geo,readvars=["pr"],
    #                                 ntimes=365,firstfile_timeinit=92,latvar="lat",lonvar="lon",subset=[10,-10,10,-10])
    #             (meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration)=runfor(wrf_data,dataareprecip=True)
    #             outprefix="wrf_dyn_12km_precip"
    #             write_data(outprefix,meanp,wet_frac,wet_duration,stddev,m90,m95,m99,maxprec99_nday,maxdry,dry_duration,reverse_index=True)
    #         else:
    #             print("No files for : wrf precip 12km")
    # 



    
    
if __name__ == '__main__':
    main()