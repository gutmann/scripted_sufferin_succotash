#!/opt/python/bin/python
###!/usr/bin/env python
import glob
import time
import cPickle
import os,re,sys
import gc
import argparse,traceback

import numpy as np

import nc_addvars
import async
import swim_io
# import mygis
# from nc import NC_writer
import date_fun
from bunch import Bunch
import fix_async_names

outputdir="output/"
geo_ref_file=None
geolat=None
geolon=None
xmin=0
xmax=None
ymin=0
ymax=None
subxmin=0
subxmax=None
subymin=0
subymax=None
geoLUT=None

# def read_geo_latlon(geo_file=None,subset=None,mask=False):
#     if not geo_file:
#         geo_file=geo_ref_file
#     if subset:
#         xmin=subset[0]
#         xmax=subset[1]
#         ymin=subset[2]
#         ymax=subset[3]
#     else:
#         # xmin=200;xmax=400;ymin=150;ymax=300;
#         xmin=0;xmax=None;ymin=0;ymax=None;
#     try:
#         lat=swim_io.read_nc(geo_file,var="lat").data
#         lon=swim_io.read_nc(geo_file,var="lon").data
#     except:
#         try:
#             lat=swim_io.read_nc(geo_file,var="latitude").data
#             lon=swim_io.read_nc(geo_file,var="longitude").data
#         except:
#             lat=swim_io.read_nc(geo_file,var="XLAT").data
#             lon=swim_io.read_nc(geo_file,var="XLONG").data
#         
#     if lon.max()>180:
#         lon=lon-360
#     if len(lon.shape)==1:
#         output=np.meshgrid(lon[xmin:xmax],lat[ymin:ymax])
#     else:
#         output=(lon[ymin:ymax,xmin:xmax],lat[ymin:ymax,xmin:xmax])
#     
#     if type(mask)==str:
#         d=swim_io.read_nc(geo_file,mask).data[0,...]
#         d=d[ymin:ymax,xmin:xmax]
#         try:
#             badpoints=np.where(d.mask)
#         except AttributeError:
#             badpoints=np.where((d<-1000) | (d>1000))
#         output[0][badpoints]=-9999
#         output[1][badpoints]=-9999
#     return output
# 
# def bilin_weights(yi,y,xi,x):
#     if (x[1]-x[0])==0:
#         x0=0.5
#     else:
#         x0=np.abs((xi-x[0])/(x[1]-x[0]))
#     x1=1-x0 # np.abs((xi-x[1])/(x[1]-x[0]))
#     if (x[3]-x[2])==0:
#         x2=0.5
#     else:
#         x2=np.abs((xi-x[2])/(x[3]-x[2]))
#     x3=1-x2 # np.abs((xi-x[3])/(x[3]-x[2]))
#     y5=y[0]*x1+y[1]*x0
#     y6=y[2]*x3+y[3]*x2
#     if (y6-y5)==0:
#         f1=0.5
#     else:
#         f1=(yi-y5)/(y6-y5)
#     f2=1-f1# (y6-yi)/(y6-y5)
#     return np.array([x1*f2,x0*f2,x3*f1,x2*f1])
# 
# def load_geoLUT(lat1,lon1,georef,subset=None,mask=False,omask=False):
#     '''Create a Geographic Look Up Table
#     
#     lat/lon inputs should be grids of latitude and longitude
#         lat1/lon1 = low resolution input grid
#         lat2/lon2 = high resolution output grid
#     '''
#     (lon2,lat2)=read_geo_latlon(georef,subset=subset,mask=mask)
#     if lon1.max()>180:
#         lon1-=360
#     N=lat2.shape
#     N1=lat1.shape
#     # output data
#     geoLUT=np.empty((N[0],N[1],4,3))
# 
#     # intermediate variables
#     x=np.zeros(4).astype('i')
#     y=np.zeros(4).astype('i')
#     winhalfsize=3 # search window for next data point
#     
#     # figure out which direction is positive in latitude and longitude (just in case)
#     if len(lon1.shape)==1:
#         lon1,lat1=np.meshgrid(lon1,lat1)
#         N1=lat1.shape
#     
#     cx=round(N1[1]/2)
#     cy=round(N1[0]/2)
#     dxinc=np.sign(lon1[cy+1,cx+1]-lon1[cy,cx]).astype('i')
#     dyinc=np.sign(lat1[cy+1,cx+1]-lat1[cy,cx]).astype('i')
# 
#     print(N1,N)
#     # loop over y
#     for i in range(N[1]):
#         # find the first positions so the rest can be done faster relative to it
#         j=0
#         # compute distances for the ENTIRE low-res grid for the first input point
#         dists=(lat1-lat2[j,i])**2 + (lon1-lon2[j,i])**2
#         # find the minimum distance position
#         (lasty,lastx)=np.unravel_index(dists.argmin(), dists.shape)
#         (prevx,prevy)=(lastx,lasty)
#         # create a window around that position to start searching as we loop through the rest of it
#         y0=max(lasty-winhalfsize,0)
#         y1=min(lasty+winhalfsize,N1[0])
#         x0=max(lastx-winhalfsize,0)
#         x1=min(lastx+winhalfsize,N1[1])
#         latwin=lat1[y0:y1,x0:x1]
#         lonwin=lon1[y0:y1,x0:x1]
#         for j in range(N[0]):
#             # if we have moved update the window position
#             if (prevx!=lastx) or (prevy!=lasty):
#                 y0=max(lasty-winhalfsize,0)
#                 y1=min(lasty+winhalfsize,N1[0])
#                 x0=max(lastx-winhalfsize,0)
#                 x1=min(lastx+winhalfsize,N1[1])
#                 latwin=lat1[y0:y1,x0:x1]
#                 lonwin=lon1[y0:y1,x0:x1]
#                 (prevx,prevy)=(lastx,lasty)
#             # compute distances just over this window
#             dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
#             # find minimum distance position
#             (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
#             # convert that position to a position in the large grid
#             lastx=newx+x0
#             lasty=newy+y0
#             # this is the first point we care about
#             x[0]=newx
#             y[0]=newy
#             # now find the position of the other three surrounding points
#             winsz=latwin.shape
#             if latwin[newy,newx]<lat2[j,i]:
#                 if lonwin[newy,newx]<lon2[j,i]:
#                     x[1]=max(min(newx+dxinc,winsz[1]-1),0)
#                     x[2]=newx
#                     x[3]=max(min(newx+dxinc,winsz[1]-1),0)
#                     y[1]=newy
#                     y[2]=max(min(newy+dyinc,winsz[0]-1),0)
#                     y[3]=max(min(newy+dyinc,winsz[0]-1),0)
#                 else:
#                     x[1]=max(min(newx-dxinc,winsz[1]-1),0)
#                     x[2]=newx
#                     x[3]=max(min(newx-dxinc,winsz[1]-1),0)
#                     y[1]=newy
#                     y[2]=max(min(newy+dyinc,winsz[0]-1),0)
#                     y[3]=max(min(newy+dyinc,winsz[0]-1),0)
#             else:
#                 if lonwin[newy,newx]<lon2[j,i]:
#                     x[1]=max(min(newx+dxinc,winsz[1]-1),0)
#                     x[2]=newx
#                     x[3]=max(min(newx+dxinc,winsz[1]-1),0)
#                     y[1]=newy
#                     y[2]=max(min(newy-dyinc,winsz[0]-1),0)
#                     y[3]=max(min(newy-dyinc,winsz[0]-1),0)
#                 else:
#                     x[1]=max(min(newx-dxinc,winsz[1]-1),0)
#                     x[2]=newx
#                     x[3]=max(min(newx-dxinc,winsz[1]-1),0)
#                     y[1]=newy
#                     y[2]=max(min(newy-dyinc,winsz[0]-1),0)
#                     y[3]=max(min(newy-dyinc,winsz[0]-1),0)
#             # finally compute the weights for each of the four surrounding points for a bilinear interpolation
#             weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x])
#             if type(mask)==np.ndarray:
#                 if np.max(x+x0)>=mask.shape[1]:
#                     tmp=np.where((x+x0)>=mask.shape[1])
#                     x[tmp]=mask.shape[1]-x0-1
#                 if np.max(y+y0)>=mask.shape[0]:
#                     tmp=np.where((y+y0)>=mask.shape[0])
#                     y[tmp]=mask.shape[0]-y0-1
#                 weights[mask[y+y0,x+x0]]=0
#             weights[weights<0]=0
#             # if there were some non-masked points, normalize by the sum of all non-masked points
#             # otherwise we will set some weights to 0 without increasing the other weights
#             if np.max(weights)>0:
#                 weights/=np.sum(weights)
#             else:
#                 if type(omask)==np.ndarray:
#                     if not omask[j,i]:
#                         sub_y0=max(lasty-(winhalfsize+3),0)
#                         sub_y1=min(lasty+(winhalfsize+3),N1[0])
#                         sub_x0=max(lastx-(winhalfsize+3),0)
#                         sub_x1=min(lastx+(winhalfsize+3),N1[1])
#                         tmp=np.where(mask[sub_y0:sub_y1,sub_x0:sub_x1]==False)
#                         if len(tmp[0])>0:
#                             latwin=lat1[sub_y0:sub_y1,sub_x0:sub_x1]
#                             lonwin=lon1[sub_y0:sub_y1,sub_x0:sub_x1]
#                             dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
# 
#                             gooddists=dists[tmp]
#                             (newy,newx)=np.where(dists==gooddists[gooddists.argmin()])
#                             x[:]=0
#                             y[:]=0
#                             x[0]=newx+sub_x0-x0
#                             y[0]=newy+sub_y0-y0
#                             weights[:]=[1,0,0,0]
#             
#             # store the current results in the output array
#             geoLUT[j,i,:,0]=y+y0
#             geoLUT[j,i,:,1]=x+x0
#             geoLUT[j,i,:,2]=weights
#     return geoLUT
#     
# def regrid(data,lat,lon,geoLUT,ymin,xmin):
#     outputdata=np.zeros((data.shape[0],geoLUT.shape[0],geoLUT.shape[1]),dtype=np.float32)
#     for i in range(4):
#         y=geoLUT[:,:,i,0].astype('i')-ymin
#         x=geoLUT[:,:,i,1].astype('i')-xmin
#         outputdata+=np.float32(data[:,y,x]*geoLUT[np.newaxis,:,:,i,2])
#     return outputdata
# 
# def convert2daily(output,days):
#     N=output.shape
#     outputdata=output.reshape((N[0]/8,8,N[1],N[2])).sum(axis=1)
#     return outputdata
# 
# def wrap_around(times,length):
#     tmp=np.where(times<0)
#     if len(tmp[0])>0:
#         times[tmp]=length+times[tmp]
#     tmp=np.where(times>=length)
#     if len(tmp[0])>0:
#         times[tmp]=times[tmp]-length
#     return times
# 
# def find_times(alltimes,month,t0,pad_length=15):
#     day0=date_fun.datearr2mjd(t0)
#     days=alltimes+day0 #convert to Modified Julian Day
#     dates=date_fun.mjd2date(days)
#     lastmonth=(month-1)%12
#     if lastmonth==0:lastmonth=12
#     nextmonth=(month+1)%12
#     if nextmonth==0:nextmonth=12
#     times_in_month=np.where((dates[:,1]==month))[0]
#     
#     # pad the data by ~two weeks on either side to make the regressions more stable
#     if pad_length>0:
#         start_padding=times_in_month[0]-(np.arange(pad_length)+1)
#         start_padding=wrap_around(start_padding,len(alltimes))
#         end_padding=times_in_month[-1]+(np.arange(pad_length)+1)
#         end_padding=wrap_around(end_padding,len(alltimes))
#         return np.hstack([start_padding,times_in_month,end_padding])
#     else:
#         return times_in_month
#         
#     
# def read_bcca_file(filename,month,geo,loadvar="prate",pad_length=15,subset=None,startdate=[1800,1,1,0,0.,0.]):
#     d=swim_io.read_nc(filename,var=loadvar)#,returnNCvar=True)
#     t=swim_io.read_nc(filename,var="time").data
#     lat=swim_io.read_nc(filename,var="lat").data
#     lon=swim_io.read_nc(filename,var="lon").data
#     subymin=int(geo[:,:,:,0].min())
#     subymax=int(geo[:,:,:,0].max()+1)
#     subxmin=int(geo[:,:,:,1].min())
#     subxmax=int(geo[:,:,:,1].max()+1)
#     
#     thismonth=find_times(t/24.0,month,np.array(startdate),pad_length=pad_length)
#     dt=(t[1]-t[0])*60*60 # t is in hours, need dt in seconds
#     # Nio doesn't let you subset with a list, only with slices?
#     output=d.data[:,subymin:subymax,subxmin:subxmax][thismonth,:,:]
#     # if using the netcdf4-python library this might be a little faster
#     # output=d.data[thismonth,subymin:subymax,subxmin:subxmax]
#     if loadvar=="prate":
#         print(dt,output.shape)
#         output*=dt #convert mm/s to mm 
#         output=convert2daily(output,t[thismonth]/24.0)
#         print(output.shape)
#     output=regrid(output,lat,lon,geo,subymin,subxmin)
#     # d.ncfile.close()
#     return output
# 
# def read_obs_file(filename,month,geo,loadvar="pr",pad_length=15,subset=None,startdate=[1940,1,1,0,0.,0.]):
#     d=swim_io.read_nc(filename,var=loadvar,returnNCvar=True)
#     t=swim_io.read_nc(filename,var="time").data
#     if subset:
#         xmin=subset[0]
#         xmax=subset[1]
#         ymin=subset[2]
#         ymax=subset[3];
#     else:
#         xmin=0;xmax=None;ymin=0;ymax=None;
#     
#     year=int(filename.split('/')[-1].split('.')[-2])
#     # if re.match(".*obs/uw.0.*",filename):
#     #     print("WARNING: correcting dates in uw-6km files, is this still necessary?")
#     #     # correct screwed up dates in the UW nc files.  
#     #     if (year%4==0):
#     #         if (year%100!=0) or (year%400==0):
#     #             t[60:]-=1
#     thismonth=find_times(t,month,np.array(startdate),pad_length=pad_length)
#     # Nio doesn't let you subset with a list, only with slices?
#     output=d.data[:,ymin:ymax,xmin:xmax][thismonth,:,:]
#     # if using the netcdf4-python library this might be a little faster
#     # output=d.data[thismonth,ymin:ymax,xmin:xmax]
#     
#     d.ncfile.close()
#     return output
#     
# 
# def read_data(files,month,filereader,geo,loadvar=None,pad_length=15,subset=None,startdate=None):
#     d1=filereader(files[0],month,geo,loadvar=loadvar,pad_length=pad_length,subset=subset,startdate=startdate)
#     output=np.empty((len(files)*(d1.shape[0]+1),d1.shape[1],d1.shape[2]))
#     output[0:d1.shape[0],:,:]=d1
#     start=d1.shape[0]
# 
#     lengths=[start]
#     for f in files[1:]:
#         if loadvar!=None:
#             d1=filereader(f,month,geo,loadvar=loadvar,pad_length=pad_length,subset=subset,startdate=startdate)
#         else:
#             d1=filereader(f,month,geo,pad_length=pad_length,subset=subset,startdate=startdate)
#         last=start+d1.shape[0]
#         output[start:last,:,:]=d1
#         start=last
#         lengths.append(start)
#         curyear=f.split('/')[-1].split('.')[1]
#         years.append(curyear)
#     output=output[:last,:,:]
#     return Bunch(data=output,dates=years,lengths=lengths)

def make_obs_times(files,data):
    """return list of MJDs corresponding to days in data"""
    dates=[]
    for f,d in zip(files,data):
        days=d.shape[0]
        yearmo=f.split(".")[-2]
        date=date_fun.date2mjd(int(yearmo[:4]),int(yearmo[4:]),1,12,0)
        dates.extend(np.arange(days)+date)
    return dates

def get_times(dates,month,year,pad):
    """docstring for get_times(dates,month,year,pad)"""
    
    goodtimes=None
    
    for n,i in enumerate(dates):
        gooddate=False
        curdate=date_fun.mjd2date(i)
        if curdate[0]==year and curdate[1]==month:
            gooddate=True
        curdate=date_fun.mjd2date(i-pad)
        if curdate[0]==year and curdate[1]==month:
            gooddate=True
        curdate=date_fun.mjd2date(i+pad)
        if curdate[0]==year and curdate[1]==month:
            gooddate=True
            
        if gooddate:
            if goodtimes==None:
                goodtimes=[]
            goodtimes.append(n)
        
    return goodtimes

def read_data(files,month,varname,pad_length,subset,startdate):
    
    data=swim_io.read_files(files,varname)
    if len(files)>10:
        dates=make_obs_times(files,data)
    else:
        times=np.concatenate(swim_io.read_files(files,"time"))
        dates=date_fun.datearr2mjd(startdate)+times
    
    data=np.concatenate(data)
    start=0
    lengths=[]
    years=[]
    for y in range(1940,2010):
        curtimes=get_times(dates,month,y,pad_length)
        if curtimes!=None:
            endpt=start+len(curtimes)
            data[start:endpt,:,:314]=data[curtimes,:,:314]
            start=endpt

            lengths.append(endpt)
            years.append(y)
    output=np.zeros((endpt,data.shape[1],314))
    output[:]=data[:endpt,:,:314]
    del data
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
        if (varname=="tasmax") or (varname=="tasmin"):
            reasonable_max=outputdata[outputdata<70].max()
            outputdata[(outputdata<10000)&(outputdata>reasonable_max)]=reasonable_max
            print(reasonable_max)
        
        swim_io.write(filename,outputdata,varname=varname)
    

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
    d=swim_io.read_nc(filename,load_var,returnNCvar=True)
    data=d.data[0,...]
    d.ncfile.close()
    if type(data.mask)==bool:
        return False
    else:
        return data.mask

def async_bcca(var=None,exp="ccsm",res="6km",forcing="BCCA",runmonth=None):
    """Perform an Asynchronous regression on bcca or NCEP data"""
    base_dir="/d2/gutmann/usbr/bcca"
    bcca_dir=base_dir+"/ccsm/"
    
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
    
    startdate=np.array([1950,1,1,0,0,0.]) 
    obs_startdate=None # np.array([1940,1,1,0,0,0.]) 
    obsextra="/*"
    obs_dir=base_dir+"/obs/"
    subset=[0,None,0,314] #Headwaters subset for 6km obs?
    
    bcca_vars=["pr"]
    obs_vars=["Prec"]

    training_search=["*19[7-9]*.nc","*200[0-5]*.nc"]
    # training_search=["*1990*.nc"]
        
    output_search="*.nc" #for all runs (should also do *.nc?)
    
    for bcca_var,obs_var in zip(bcca_vars,obs_vars):
        bcca_var_dir=obs_var+"/"
        bccafiles=[]
        load_bcca_var="pr"
        
        for cur_search in training_search:
            bccafiles.extend(glob.glob(bcca_dir+cur_search))
        bccafiles.sort()
        bccafiles=list(np.unique(bccafiles))
        outputfiles=glob.glob(bcca_dir+output_search)
        outputfiles.sort()
        
        obs_files=[]
        for cur_search in training_search:
            obs_files.extend(glob.glob(obs_dir+cur_search))
        obs_files.sort()
        
        # geo_ref_file=obs_files[0]
        # (lon,lat)=read_geo_latlon(bccafiles[0])#,mask=load_bcca_var)

        # print("Creating Geographic Look Up Table")
        # t1=time.time()
        # mask=read_mask(bccafiles[0],load_bcca_var)
        # omask=read_mask(obs_files[0],obs_var)
        # geo=load_geoLUT(lat,lon,geo_ref_file,subset=subset,mask=mask,omask=omask)
        # print("   "+str((time.time()-t1)/60)+" minutes")
        
        for month in range(startmonth,endmonth):
            print("Month: "+str(month))
            geo=None
            print("  Loading "+forcing+" data")
            t1=time.time()
            bcca=read_data(bccafiles,month,load_bcca_var,pad_length=15,subset=subset,startdate=startdate)
            print("   "+str((time.time()-t1)/60)+" minutes")
            sys.stdout.flush()
        
            print("  Loading obs data")
            t1=time.time()
            obs=read_data(obs_files,month,obs_var,pad_length=15,subset=subset,startdate=obs_startdate)
            print("   "+str((time.time()-t1)/60)+" minutes")
        
            print("  Computing async regression")
            t1=time.time()
            # just in case one data set is longer than the other
            print("Obs data length:",obs.data.shape[0])
            print("BCCA data length:",bcca.data.shape[0])
            tmin=min(obs.data.shape[0],bcca.data.shape[0])
            obsmax=obs.data.max()
            regression=async.develop_async(hires=obs.data[:tmin,...],lowres=bcca.data[:tmin,...],
                    isPrecip=(bcca_var=="pr"),verbose=True,even_xy=False)
            print("   "+str((time.time()-t1)/60)+" minutes")
         
            print("  Applying async regression")
            t1=time.time() 
            bcca=read_data(outputfiles,month,load_bcca_var,pad_length=0,subset=subset,startdate=startdate)
            output=async.apply_async(bcca.data,regression,vmax=obsmax*1.2,isPrecip=(bcca_var=="pr"),verbose=True)
            # if type(omask)==np.ndarray:
            #     for thistime in range(output.shape[0]):
            #         output[thistime,:,:][omask]=-9999
        
            prefix='0' if month<10 else ''
            picklefile=open("_".join(["SAR",forcing,exp,res,bcca_var,prefix+str(month)])+".pickle",'wb')
            pickle=cPickle.Pickler(picklefile,2)
            pickle.dump(regression)
            picklefile.close()
        
            print("  Writing data")
            resextra=""
            write_data(output,month,bcca.dates,bcca.lengths,obs_var,res+resextra,output_dir=output_dir)
            
            # trying to get rid of a weird memory leak
            # cleanup large variables at the end of each month of processing. 
            del bcca
            del output
            del obs
            gc.collect()
        
        if runmonth==None:
            # add the lat and lon variables to the orignal files (should also add time?)
            nc_addvars.main(output_dir+"BCSAR*"+res+"*.nc",copy_from_file=bccafiles[0],
                            ranges=[[subset[0],subset[1]],[subset[2],subset[3]]],
                            vars2copy=["longitude","latitude"])
            # next delete the old files
            old_files=glob.glob(output_dir+"BCSAR*"+res+"*.nc")
            old_files.sort()
            for f in old_files:
                print(f)
                os.remove(f)
            # then rename "added" files back to the original filenames
            new_files=glob.glob(output_dir+"added*"+res+"*.nc")
            new_files.sort()
            for o,n in zip(old_files,new_files):
                print(n,o)
                os.rename(n,o)
            
        
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Compute a Statistical Asynchronous Regression (with hard coded files). ')
        parser.add_argument('var',action='store',nargs="?",help="Chose a variable [tasmax,tasmin,pr]",default="pr")
        parser.add_argument('exp',action='store',nargs="?",help="Chose an experiment [e0,e1,conus,wet,dry,hot,cold]",default="e0")
        parser.add_argument('res',action='store',nargs="?",help="Chose a resolution [4km, 6km, 12km]",default="6km")
        parser.add_argument('forcing',action='store',nargs="?",help="Chose a forcing model [NCEP, bcca]",default="bcca")
        parser.add_argument('month',action='store',nargs="?",help="Only run this month",default=None)
        parser.add_argument('-v', '--version',action='version',
                version='async_bcca.py 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
    
        exit_code = async_bcca(args.var,args.exp,args.res,args.forcing,args.month)
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
    