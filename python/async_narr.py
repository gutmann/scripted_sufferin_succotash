#!/usr/bin/env python
import glob
import time
import cPickle
import os,re,sys
import gc
import argparse,traceback

import numpy as np

import async
import mygis
import date_fun
from bunch import Bunch

# from flushprint import Flushfile
# sys.stdout=Flushfile(sys.stdout)


# outputdir="output/"
# geo_ref_file=None
# xmin=0
# xmax=None
# ymin=0
# ymax=None
# subxmin=0
# subxmax=None
# subymin=0
# subymax=None
# geoLUT=None

# CCSM specific offsets into a single huge file
startyears =Bunch(train=1979,apply=1979)
endyears   =Bunch(train=1999,apply=2060)

ccsm_start_year=1900
# ccsm_start_year=1960
startpoints=Bunch(train=(startyears["train"]-ccsm_start_year)*365,
                  apply=(startyears["apply"]-ccsm_start_year)*365)
endpoints=Bunch(train=(endyears["train"]-ccsm_start_year+1)*365,
                  apply=(endyears["apply"]-ccsm_start_year+1)*365)

days_per_month=[31,28,31,30,31,30,31,31,30,31,30,31]

def read_geo_latlon(geo_file=None,subset=None,mask=False):
    if not geo_file:
        geo_file=geo_ref_file
    if subset:
        xmin=subset[0]
        xmax=subset[1]
        ymin=subset[2]
        ymax=subset[3]
    else:
        # xmin=200;xmax=400;ymin=150;ymax=300;
        xmin=0;xmax=None;ymin=0;ymax=None;
    try:
        lat=mygis.read_nc(geo_file,var="lat").data
        lon=mygis.read_nc(geo_file,var="lon").data
    except:
        try:
            lat=mygis.read_nc(geo_file,var="latitude").data
            lon=mygis.read_nc(geo_file,var="longitude").data
        except:
            lat=mygis.read_nc(geo_file,var="XLAT").data
            lon=mygis.read_nc(geo_file,var="XLONG").data
        
    if lon.max()>180:
        lon=lon-360
    if len(lon.shape)==1:
        output=np.meshgrid(lon[xmin:xmax],lat[ymin:ymax])
    else:
        output=(lon[ymin:ymax,xmin:xmax],lat[ymin:ymax,xmin:xmax])
    
    if type(mask)==str:
        d=mygis.read_nc(geo_file,mask).data[0,...]
        d=d[ymin:ymax,xmin:xmax]
        try:
            badpoints=np.where(d.mask)
        except AttributeError:
            badpoints=np.where((d<-1000) | (d>1000))
        output[0][badpoints]=-9999
        output[1][badpoints]=-9999
    return output

def bilin_weights(yi,y,xi,x):
    if (x[1]-x[0])==0:
        x0=0.5
    else:
        x0=np.abs((xi-x[0])/(x[1]-x[0]))
    x1=1-x0 # np.abs((xi-x[1])/(x[1]-x[0]))
    if (x[3]-x[2])==0:
        x2=0.5
    else:
        x2=np.abs((xi-x[2])/(x[3]-x[2]))
    x3=1-x2 # np.abs((xi-x[3])/(x[3]-x[2]))
    y5=y[0]*x1+y[1]*x0
    y6=y[2]*x3+y[3]*x2
    if (y6-y5)==0:
        f1=0.5
    else:
        f1=(yi-y5)/(y6-y5)
    f2=1-f1# (y6-yi)/(y6-y5)
    return np.array([x1*f2,x0*f2,x3*f1,x2*f1])

def load_geoLUT(lat1,lon1,georef,subset=None,mask=False,omask=False):
    '''Create a Geographic Look Up Table
    
    lat/lon inputs should be grids of latitude and longitude
        lat1/lon1 = low resolution input grid
        lat2/lon2 = high resolution output grid
    '''
    (lon2,lat2)=read_geo_latlon(georef,subset=subset,mask=mask)
    if lon1.max()>180:
        lon1-=360
    N=lat2.shape
    N1=lat1.shape
    # output data
    geoLUT=np.empty((N[0],N[1],4,3))

    # intermediate variables
    x=np.zeros(4).astype('i')
    y=np.zeros(4).astype('i')
    winhalfsize=3 # search window for next data point
    
    # figure out which direction is positive in latitude and longitude (just in case)
    if len(lon1.shape)==1:
        lon1,lat1=np.meshgrid(lon1,lat1)
        N1=lat1.shape
    
    cx=round(N1[1]/2)
    cy=round(N1[0]/2)
    dxinc=np.sign(lon1[cy+1,cx+1]-lon1[cy,cx]).astype('i')
    dyinc=np.sign(lat1[cy+1,cx+1]-lat1[cy,cx]).astype('i')

    print(N1,N)
    # loop over y
    for i in range(N[1]):
        # find the first positions so the rest can be done faster relative to it
        j=0
        # compute distances for the ENTIRE low-res grid for the first input point
        dists=(lat1-lat2[j,i])**2 + (lon1-lon2[j,i])**2
        # find the minimum distance position
        (lasty,lastx)=np.unravel_index(dists.argmin(), dists.shape)
        (prevx,prevy)=(lastx,lasty)
        # create a window around that position to start searching as we loop through the rest of it
        y0=max(lasty-winhalfsize,0)
        y1=min(lasty+winhalfsize,N1[0])
        x0=max(lastx-winhalfsize,0)
        x1=min(lastx+winhalfsize,N1[1])
        latwin=lat1[y0:y1,x0:x1]
        lonwin=lon1[y0:y1,x0:x1]
        for j in range(N[0]):
            # if we have moved update the window position
            if (prevx!=lastx) or (prevy!=lasty):
                y0=max(lasty-winhalfsize,0)
                y1=min(lasty+winhalfsize,N1[0])
                x0=max(lastx-winhalfsize,0)
                x1=min(lastx+winhalfsize,N1[1])
                latwin=lat1[y0:y1,x0:x1]
                lonwin=lon1[y0:y1,x0:x1]
                (prevx,prevy)=(lastx,lasty)
            # compute distances just over this window
            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
            # find minimum distance position
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            # convert that position to a position in the large grid
            lastx=newx+x0
            lasty=newy+y0
            # this is the first point we care about
            x[0]=newx
            y[0]=newy
            # now find the position of the other three surrounding points
            winsz=latwin.shape
            if latwin[newy,newx]<lat2[j,i]:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=max(min(newx+dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx+dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy+dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy+dyinc,winsz[0]-1),0)
                else:
                    x[1]=max(min(newx-dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx-dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy+dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy+dyinc,winsz[0]-1),0)
            else:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=max(min(newx+dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx+dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy-dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy-dyinc,winsz[0]-1),0)
                else:
                    x[1]=max(min(newx-dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx-dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy-dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy-dyinc,winsz[0]-1),0)
            # finally compute the weights for each of the four surrounding points for a bilinear interpolation
            weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x])
            if type(mask)==np.ndarray:
                if np.max(x+x0)>=mask.shape[1]:
                    tmp=np.where((x+x0)>=mask.shape[1])
                    x[tmp]=mask.shape[1]-x0-1
                if np.max(y+y0)>=mask.shape[0]:
                    tmp=np.where((y+y0)>=mask.shape[0])
                    y[tmp]=mask.shape[0]-y0-1
                weights[mask[y+y0,x+x0]]=0
            weights[weights<0]=0
            # if there were some non-masked points, normalize by the sum of all non-masked points
            # otherwise we will set some weights to 0 without increasing the other weights
            if np.max(weights)>0:
                weights/=np.sum(weights)
            else:
                if type(omask)==np.ndarray:
                    if not omask[j,i]:
                        sub_y0=max(lasty-(winhalfsize+3),0)
                        sub_y1=min(lasty+(winhalfsize+3),N1[0])
                        sub_x0=max(lastx-(winhalfsize+3),0)
                        sub_x1=min(lastx+(winhalfsize+3),N1[1])
                        tmp=np.where(mask[sub_y0:sub_y1,sub_x0:sub_x1]==False)
                        if len(tmp[0])>0:
                            latwin=lat1[sub_y0:sub_y1,sub_x0:sub_x1]
                            lonwin=lon1[sub_y0:sub_y1,sub_x0:sub_x1]
                            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2

                            gooddists=dists[tmp]
                            (newy,newx)=np.where(dists==gooddists[gooddists.argmin()])
                            x[:]=0
                            y[:]=0
                            x[0]=newx+sub_x0-x0
                            y[0]=newy+sub_y0-y0
                            weights[:]=[1,0,0,0]
            
            # store the current results in the output array
            geoLUT[j,i,:,0]=y+y0
            geoLUT[j,i,:,1]=x+x0
            geoLUT[j,i,:,2]=weights
    return geoLUT
    
def regrid(data,lat,lon,geoLUT,ymin,xmin):
    outputdata=np.zeros((data.shape[0],geoLUT.shape[0],geoLUT.shape[1]),dtype=np.float32)
    for i in range(4):
        y=geoLUT[:,:,i,0].astype('i')-ymin
        x=geoLUT[:,:,i,1].astype('i')-xmin
        outputdata+=np.float32(data[:,y,x]*geoLUT[np.newaxis,:,:,i,2])
    return outputdata

def convert2daily(output,days):
    N=output.shape
    outputdata=output.reshape((N[0]/8,8,N[1],N[2])).sum(axis=1)
    return outputdata

def wrap_around(times,length):
    # checks for areas where times is outside of the bounds 0,length
    # wraps times such that times >length start counting from 0 again
    #  times less than 0 start counting down from length
    tmp=np.where(times<0)
    if len(tmp[0])>0:
        times[tmp]=length+times[tmp]
    tmp=np.where(times>=length)
    if len(tmp[0])>0:
        times[tmp]=times[tmp]-length
    return times

def find_times(alltimes,month,t0,pad_length=15):
    day0=date_fun.datearr2mjd(t0)
    days=alltimes+day0 #convert to Modified Julian Day
    dates=date_fun.mjd2date(days)
    
    times_in_month=np.where((dates[:,1]==month))[0]
    
    # pad the data by ~two weeks on either side to make the regressions more stable
    if pad_length>0:
        curtime=times_in_month[0]
        output_times=times_in_month.copy()
        additions=0
        for i in range(1,len(times_in_month)-pad_length):
            if (times_in_month[i]-curtime)>min(pad_length,359):
                # we just jumped years
                print(i,additions,output_times.shape,pad_length)
                new_output_times=np.zeros(output_times.shape[0]+pad_length*2,dtype=np.int)
                new_output_times[:i+additions]=output_times[:i+additions]
                new_output_times[i+additions:i+additions+pad_length]=output_times[i+additions-1]+np.arange(1,pad_length+1)
                new_output_times[i+additions+pad_length:i+additions+pad_length*2]=output_times[i+additions]+np.arange(-pad_length,0,1)
                new_output_times[i+additions+pad_length*2:]=output_times[i+additions:]
                output_times=new_output_times
                additions+=pad_length*2
            curtime=times_in_month[i]
        times_in_month=output_times
                
        start_padding=times_in_month[0]-(np.arange(pad_length)+1)
        start_padding=wrap_around(start_padding,len(alltimes))
        end_padding=times_in_month[-1]+(np.arange(pad_length)+1)
        end_padding=wrap_around(end_padding,len(alltimes))
        return np.hstack([start_padding,times_in_month,end_padding])
        
    else:
        return times_in_month


def add_trend_back(data,trend_start_point,trend,points_per_year=365):
    n=data.shape[0]-trend_start_point
    if n<(points_per_year*2):
        # we have less than 2 years of data, don't bother
        return
    x=np.arange(n)
    # plt.plot(data[trend_start_point:,150:189,100:230].mean(axis=1).mean(axis=1),label="AR-pre-addition")
    data[trend_start_point:,:,:]+=x[:,np.newaxis,np.newaxis]*trend.gain[np.newaxis,:,:]
    # plt.plot(data[trend_start_point:,150:189,100:230].mean(axis=1).mean(axis=1),label="AR-post-addition")
    # plt.legend()
    # plt.savefig("test.png")

def remove_trend_in_data(data,trend_start_point,points_per_year=365):
    n=data.shape[0]-trend_start_point
    if n<(points_per_year*2):
        # we have less than 2 years of data, don't bother
        return
    start_means=data[trend_start_point:trend_start_point+365,:,:].mean(axis=0)
    # plt.clf()
    # plt.plot(data[trend_start_point:,150:189,100:230].mean(axis=1).mean(axis=1),label="CCSM-pre")
    data[trend_start_point:,:,:]-=start_means[np.newaxis,:,:]
    mid_point_means=data[trend_start_point:,:,:].mean(axis=0)/(n/2.0)
    x=np.arange(n)
    data[trend_start_point:,:,:]-=x[:,np.newaxis,np.newaxis]*mid_point_means[np.newaxis,:,:]
    data[trend_start_point:,:,:]+=start_means[np.newaxis,:,:]
    # plt.plot(data[trend_start_point:,150:189,100:230].mean(axis=1).mean(axis=1),label="CCSM-post")
    return Bunch(offset=start_means,gain=mid_point_means)
        
    
def read_narr_file(filename,month,geo,loadvar="prate",pad_length=15,subset=None,startdate=[1800,1,1,0,0.,0.]):
    d=mygis.read_nc(filename,var=loadvar)#,returnNCvar=True)
    t=mygis.read_nc(filename,var="time").data
    lat=mygis.read_nc(filename,var="lat").data
    lon=mygis.read_nc(filename,var="lon").data
    subymin=int(geo[:,:,:,0].min())
    subymax=int(geo[:,:,:,0].max()+1)
    subxmin=int(geo[:,:,:,1].min())
    subxmax=int(geo[:,:,:,1].max()+1)
    
    thismonth=find_times(t/24.0,month,np.array(startdate),pad_length=pad_length)
    dt=(t[1]-t[0])*60*60 # t is in hours, need dt in seconds
    # Nio doesn't let you subset with a list, only with slices?
    output=d.data[:,subymin:subymax,subxmin:subxmax][thismonth,:,:]
    # if using the netcdf4-python library this might be a little faster
    # output=d.data[thismonth,subymin:subymax,subxmin:subxmax]
    if loadvar=="prate":
        print(dt,output.shape)
        output*=dt #convert mm/s to mm 
        output=convert2daily(output,t[thismonth]/24.0)
        print(output.shape)
    output=regrid(output,lat,lon,geo,subymin,subxmin)
    # d.ncfile.close()
    return output

def read_ccsm_file(filename,month,geo,loadvar="pr",pad_length=15,subset=None,startdate=[0,1,1,0,0.,0.],period="train"):
    t=mygis.read_nc(filename,var="time").data
    lat=mygis.read_nc(filename,var="lat").data
    lon=mygis.read_nc(filename,var="lon").data
    d=mygis.read_nc(filename,var=loadvar,returnNCvar=True)

    startpoint=startpoints[period]
    endpoint=endpoints[period]
    # CCSM precip and tas actually start in 1960, not 1900, so subtract 60yrs (in days)
    if ((loadvar=="pr") or (loadvar=="tas")) and (ccsm_start_year==1900):
        print("short circuiting...")
        startpoint-=60*365
        endpoint-=60*365
    subymin=int(geo[:,:,:,0].min())
    subymax=int(geo[:,:,:,0].max()+1)
    subxmin=int(geo[:,:,:,1].min())
    subxmax=int(geo[:,:,:,1].max()+1)
    
    # because times are from a 365 day year (noleap) we have to trick find_times by using the day of year
    # this "pretends" they are all in the same year, but all find_times uses is the month so as long as
    # startdate is a non-leap year this should work
    doy=t%365
    fake_start_date=np.array([1999,1,1,0,0,0])
    thismonth=find_times(doy[startpoint:endpoint],month,fake_start_date,pad_length=pad_length)
    dt=(t[1]-t[0])*60*60*24 # t is in days, need dt in seconds (= 86400)
    
    # Nio doesn't let you subset with a list, only with slices?
    output=d.data[startpoint:endpoint,subymin:subymax,subxmin:subxmax][thismonth,:,:]
    # if using the netcdf4-python library this might be a little faster
    # output=d.data[thismonth,subymin:subymax,subxmin:subxmax]
    if loadvar=="pr":
        output*=dt #convert mm/s to mm
    output=regrid(output,lat,lon,geo,subymin,subxmin)
    d.ncfile.close()
    years=np.arange(startyears[period],endyears[period]+1)
    dayspermonth=(output.shape[0]/len(years))-pad_length*2
    print("Days in the month:",dayspermonth)
    print("Days found:",output.shape[0])
    print("Number of years:",len(years))
    return Bunch(data=output,dates=years,lengths=[dayspermonth*i+dayspermonth for i in range(len(years))])



def read_obs_file(filename,month,geo,loadvar="pr",pad_length=15,subset=None,startdate=[1940,1,1,0,0.,0.]):
    d=mygis.read_nc(filename,var=loadvar,returnNCvar=True)
    t=mygis.read_nc(filename,var="time").data
    if subset:
        xmin=subset[0]
        xmax=subset[1]
        ymin=subset[2]
        ymax=subset[3];
    else:
        xmin=0;xmax=None;ymin=0;ymax=None;
    
    year=int(filename.split('/')[-1].split('.')[-2])
    thismonth=find_times(t,month,np.array(startdate),pad_length=pad_length)
    # Nio doesn't let you subset with a list, only with slices?
    output=d.data[:,ymin:ymax,xmin:xmax][thismonth,:,:]
    # if using the netcdf4-python library this might be a little faster
    # output=d.data[thismonth,ymin:ymax,xmin:xmax]
    
    d.ncfile.close()
    return output
    

def read_data(files,month,filereader,geo,loadvar=None,pad_length=15,subset=None,startdate=None):
    
    d1=filereader(files[0],month,geo,loadvar=loadvar,pad_length=pad_length,subset=subset,startdate=startdate)
    output=np.empty((len(files)*(d1.shape[0]+1),d1.shape[1],d1.shape[2]))
    output[0:d1.shape[0],:,:]=d1
    start=d1.shape[0]
    curyear=files[0].split('/')[-1].split('.')[1]
    if curyear=="2m":
        curyear=files[0].split('/')[-1].split('.')[2]
    if curyear=="sfc":
        curyear=files[0].split('/')[-1].split('.')[-3]
    if curyear=="gauss":
        curyear=files[0].split('/')[-1].split('.')[-2]
    if curyear=="sig995":
        curyear=files[0].split('/')[-1].split('.')[-2]
    years=[curyear]
    lengths=[start]
    for f in files[1:]:
        if loadvar!=None:
            d1=filereader(f,month,geo,loadvar=loadvar,pad_length=pad_length,subset=subset,startdate=startdate)
        else:
            d1=filereader(f,month,geo,pad_length=pad_length,subset=subset,startdate=startdate)
        last=start+d1.shape[0]
        output[start:last,:,:]=d1
        start=last
        lengths.append(start)
        curyear=f.split('/')[-1].split('.')[1]
        if curyear=="2m":
            curyear=f.split('/')[-1].split('.')[2]
        if curyear=="sfc":
            curyear=f.split('/')[-1].split('.')[-3]
        if curyear=="gauss":
            curyear=f.split('/')[-1].split('.')[-2]
        if curyear=="sig995":
            curyear=f.split('/')[-1].split('.')[-2]
        years.append(curyear)
    output=output[:last,:,:]
    return Bunch(data=output,dates=years,lengths=lengths)


def write_data(data,month,years,endpts,varname,res,output_dir="async_output/"):
    n=len(years)
    if not glob.glob(output_dir):
        try:
            os.mkdir(output_dir)
        except Exception as e:
            print(e)
    startpt=0
    if month<10:
        month_prefix="0"
    else: 
        month_prefix=""
    for (year,stoppt) in zip(years,endpts):
        if (stoppt-startpt)>35:
            outputdata=data[startpt+14:stoppt-14,:,:]
        else:
            outputdata=data[startpt:stoppt,:,:]
        startpt=stoppt
        filename=output_dir+"BCSAR_"+varname+"_"+res+"_"+str(year)+"_"+month_prefix+str(month)# +".nc"
        # if (varname=="tasmax") or (varname=="tasmin") or (varname=="tas"):
        #     reasonable_max=outputdata[outputdata<70].max()
        #     outputdata[(outputdata<10000)&(outputdata>reasonable_max)]=reasonable_max
        #     print(reasonable_max)
        
        mygis.write(filename,outputdata,varname=varname)
    

def mk_all_dirs(dirname):
    subdirs=dirname.split("/")
    for curdir in subdirs:
        if os.path.isfile(curdir):
            os.rename(curdir,os.tempnam("./",curdir+"_"))
        if not os.path.isdir(curdir):
            try:
                os.mkdir(curdir)
            except Exception as e:
                print(e)
        os.chdir(curdir)
    for i in range(len(subdirs)):
        os.chdir("../")

def read_mask(filename,load_var):
    d=mygis.read_nc(filename,load_var,returnNCvar=True)
    data=d.data[0,...]
    d.ncfile.close()
    try:
        if type(data.mask)==bool:
            return False
        else:
            return data.mask
    except:
        return False

def async_narr(var=None,exp="e0",res="12km",forcing="NCEP",runmonth=None,remove_trend=False):
    """Perform an Asynchronous regression on NARR or NCEP data"""
    base_dir="/d2/gutmann/usbr/stat_data/"
    # note you also have to add "nldas" between the / and the * after obsextra
    narr_dir=base_dir+"forcing/narr/"
    narr_extra=""
    ccsm_dir="ccsm/"
    ccsm_extra=""
    ncep_dir=base_dir+"forcing/ncep/"
    ncep_extra="*gauss*"

    if runmonth==None:
        startmonth=1
        endmonth=13
        runmonth=''
    else:
        startmonth=int(runmonth)
        endmonth=int(runmonth)+1
        

    if res[-2:]!="km":res+="km"
    output_dir="/".join(["SAR3"+exp,forcing.lower(),var])
    mk_all_dirs(output_dir)
    output_dir+="/"
    if forcing=="NCEP":
        extra=ncep_extra
        narr_dir=ncep_dir
        startdate=np.array([1,1,1,0,0.,0.])
    elif forcing=="CCSM":
        extra=ccsm_extra
        narr_dir=ccsm_dir
        startdate=np.array([0,1,1,0,0.,0.])
    else:
        extra=narr_extra
        startdate=np.array([1800,1,1,0,0.,0.])
    
    obs_startdate=np.array([1940,1,1,0,0,0.])
    obsextra="/*"
    # subset=[xmin,xmax,ymin,ymax]
    if res=="12km":
        obs_dir=base_dir+"DAILY/obs/maurer.125/"
        # obsextra="/nldas*"
        subset=[75,200,60,152] #Headwaters subset for 12km obs
    elif res=="6km":
        obs_dir=base_dir+"DAILY/obs/uw.0625/"
        # obsextra="/*"
        subset=[120,400,140,300] #Headwaters subset for 6km obs?
    elif res=="4km":
        obs_dir=base_dir+"DAILY/obs/uw.4km/"
        # obsextra="/*"
        subset=[220,600,220,500] #Headwaters subset for 4km obs
        
    narr_vars=["tmin","tmax","prate"]
    obs_vars=["tasmin","tasmax","pr"]
    # narr_vars=["prate"]
    # obs_vars=["pr"]
    if var!=None:
        obs_vars=[var]
        if var=="tasmin":narr_vars=["tmin"]
        if var=="tasmax":narr_vars=["tmax"]
        if var=="tas":narr_vars=["tas"]
        if var=="pr":narr_vars=["prate"]
    hot_search=["*198[8,6,0,9,1]*.nc","*199[5,0,6,4,9]*.nc"]
    cold_search=["*1979*.nc","*198[4,2,3,5]*.nc","*199[3,7,1,2,8]*.nc"]
    wet_search=["*1979*.nc","*198[9,7,8,0]*.nc","*199[4,2,0,9,6]*.nc"]
    dry_search=["*198[5,6,4,2,3]*.nc","*199[1,3,7,8,5]*.nc"]
    e0_search=["*19*.nc"] #for e0 and "normal" runs
    e1_search=["*20*.nc"] #for e1 runs for narr Tair should subset out 2006.5-2008.0
    if exp=="e0":training_search=e0_search
    if exp=="e1":training_search=e1_search
    if exp=="wet":training_search=wet_search
    if exp=="dry":training_search=dry_search
    if exp=="hot":training_search=hot_search
    if exp=="cold":training_search=cold_search
    if (exp=="conus") or (exp=="pgw"):
        training_search=e0_search
        subset=[0,None,0,None] #no subset, =CONUS
        
    # output_search="*20*.nc" #for 20C
    output_search="*.nc" #for all runs (should also do *.nc?)
    # subset=[lon0,lon1,lat0,lat1]
    # subset=[240,245,250,255] #small/fast testing subset on 4km or 6km
    # subset=[140,145,150,155] #small/fast testing subset on 12km or 6km
    # training_search=["*199*"] #faster than e0
    print(narr_vars, obs_vars,exp,training_search)
    for narr_var,obs_var in zip(narr_vars,obs_vars):
        narr_var_dir=obs_var+"/"
        narrfiles=[]
        if forcing=="CCSM":
            if narr_var=="tmin":
                load_narr_var="TREFMN"
            if narr_var=="tmax":
                load_narr_var="TREFMX"
        else:
            if narr_var=="tmin":
                load_narr_var="tasmin"
            if narr_var=="tmax":
                load_narr_var="tasmax"

        if narr_var=="tas":
            load_narr_var="tas"
        if narr_var=="prate":
            load_narr_var="pr"
        for cur_search in training_search:
            print(narr_dir+narr_var_dir+narr_var+extra+cur_search)
            narrfiles.extend(glob.glob(narr_dir+narr_var_dir+narr_var+extra+cur_search))
        narrfiles.sort()
        if exp=="pgw":
            outputfiles=glob.glob("pgw/ncep/"+narr_var_dir+narr_var+extra+output_search)
        else:
            outputfiles=glob.glob(narr_dir+narr_var_dir+narr_var+extra+output_search)

        outputfiles.sort()
        obs_files=[]
        for cur_search in training_search:
            print(obs_dir+obs_var+obsextra+obs_var+cur_search)
            obs_files.extend(glob.glob(obs_dir+obs_var+obsextra+obs_var+cur_search))
        obs_files.sort()
        geo_ref_file=obs_files[0]
        print(narrfiles[0])
        (lon,lat)=read_geo_latlon(narrfiles[0])#,mask=load_narr_var)

        print("Creating Geographic Look Up Table")
        t1=time.time()
        mask=read_mask(narrfiles[0],load_narr_var)
        omask=read_mask(obs_files[0],obs_var)
        geo=load_geoLUT(lat,lon,geo_ref_file,subset=subset,mask=mask,omask=omask)
        print("   "+str((time.time()-t1)/60)+" minutes")
        
        for month in range(startmonth,endmonth):
            print("Month: "+str(month))
        
            print("  Loading "+forcing+" data")
            t1=time.time()
            print(narrfiles[0])
            if forcing=="CCSM":
                narr=read_ccsm_file(narrfiles[0],month,geo,load_narr_var,pad_length=15,startdate=startdate,period="train")
            else:
                narr=read_data(narrfiles,month,read_narr_file,geo,load_narr_var,pad_length=15,startdate=startdate)
            print("   "+str((time.time()-t1)/60)+" minutes")
            sys.stdout.flush()
        
            print("  Loading obs data")
            t1=time.time()
            print(obs_files[0])
            obs=read_data(obs_files,month,read_obs_file,geo,obs_var,pad_length=15,subset=subset,startdate=obs_startdate)
            print("   "+str((time.time()-t1)/60)+" minutes")
        
            print("  Computing async regression")
            t1=time.time()
            # just in case one data set is longer than the other
            tmin=min(obs.data.shape[0],narr.data.shape[0])
            obsmax=obs.data[obs.data<1e5].max()
            regression=async.develop_async(hires=obs.data[:tmin,...],lowres=narr.data[:tmin,...],
                    isPrecip=(narr_var=="prate"),verbose=True,even_xy=False)
            print("   "+str((time.time()-t1)/60)+" minutes")
         
            print("  Applying async regression")
            print("      Loading "+forcing+" data")
            t1=time.time()
            if forcing=="CCSM":
                narr=read_ccsm_file(narrfiles[0],month,geo,load_narr_var,pad_length=0,startdate=startdate,period="apply")
            else:
                narr=read_data(outputfiles,month,read_narr_file,geo,load_narr_var,pad_length=0,subset=subset,startdate=startdate)
                
            print("      "+str((time.time()-t1)/60)+" minutes")
            print("      Applying")
            t1=time.time()
            if remove_trend:
                if load_narr_var=="pr":
                    raise Exception("Removal of trend for precip not allowed")
                print("      Removing trend")
                points_per_year=days_per_month[month-1]
                if forcing=="CCSM":
                    trend_start_point=(points_per_year*(2000-ccsm_start_year)-((startpoints["apply"]/365)*points_per_year))
                else:
                    trend_start_point=0
                trend=remove_trend_in_data(narr.data,trend_start_point,points_per_year=points_per_year)
            output=async.apply_async(narr.data,regression,vmax=obsmax*1.2,isPrecip=(narr_var=="prate"),verbose=True)
            
            if remove_trend:
                add_trend_back(output,trend_start_point,trend,points_per_year=points_per_year)
                
            
            if (type(omask)==np.ndarray) and (omask.shape[0]==output.shape[1]) and (omask.shape[1]==output.shape[2]):
                for thistime in range(output.shape[0]):
                    output[thistime,:,:][omask]=-9999
        
            prefix='0' if month<10 else ''
            picklefile=open("_".join(["SAR",forcing,exp,res,narr_var,prefix+str(month)])+".pickle",'wb')
            pickle=cPickle.Pickler(picklefile,2)
            pickle.dump(regression)
            picklefile.close()
        
            print("  Writing data")
            resextra=""
            if forcing=="NCEP":
                resextra="_gauss"
            write_data(output,month,narr.dates,narr.lengths,obs_var,res+resextra,output_dir=output_dir)
                        
        
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Compute a Statistical Asynchronous Regression (with hard coded files). ')
        parser.add_argument('var',action='store',nargs="?",help="Chose a variable [tasmax,tasmin,tas,pr]",default="pr")
        parser.add_argument('exp',action='store',nargs="?",help="Chose an experiment [e0,e1,conus,pgw,wet,dry,hot,cold]",default="e0")
        parser.add_argument('res',action='store',nargs="?",help="Chose a resolution [4km, 6km, 12km]",default="12km")
        parser.add_argument('forcing',action='store',nargs="?",help="Chose a forcing model [NCEP, NARR,CCSM]",default="NCEP")
        parser.add_argument('month',action='store',nargs="?",help="Only run this month",default=None)
        parser.add_argument('--remove_trend',dest='remove_trend',action='store_true',
                help="Remove the trend over the time series",default=False)
        parser.add_argument('-v', '--version',action='version',
                version='async_narr.py 1.2')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
    
        exit_code = async_narr(args.var,args.exp,args.res,args.forcing,args.month,args.remove_trend)
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
    