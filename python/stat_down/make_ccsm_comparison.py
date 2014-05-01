#!/usr/bin/env python
import glob,sys,datetime
import numpy as np
import matplotlib.pyplot as plt

import mygis as io
from bunch import Bunch
from bilin_regrid import load_geoLUT
from stat_down import map_vis

from flushprint import Flushfile
sys.stdout=Flushfile(sys.stdout)

run="tas"

geofile="4km_HW_metadata.nc"
ccsm_var_to_read="PRCP"
ccsm_cur="ccsm-delta/b30.030e.cam2.h0.avgPRCP.1995-01_cat_2005-12.nc"
ccsm_fut="ccsm-delta/b30.042e.cam2.h0.avgPRCP.2045-01_cat_2055-12.nc"

if run=="pr":
    esgfile="ccsm/pr/prate.run5.complete.1960-2080.nc"
    esgvarname="pr"
else:
    esgfile="ccsm/tas/tas.run5.complete.1960-2080.nc"
    esgvarname="tas"

wrf_cur=dict(pr="wrf/ctrl/NAR*",tas="wrf/ctrl/WRF_T2m-ctl.nc")
wrf_fut=dict(pr="wrf/pgw/PGW*",tas="wrf/pgw/WRF_T2m-pgw.nc")


ar_current_list=dict(pr=["sar/pr/*199[5-9]","sar/pr/*200[0-5]"],tas=["sar/tas/*199[5-9]","sar/tas/*200[0-5]"])
ar_future_list=dict(pr=["sar/pr/*204[5-9]","sar/pr/*205[0-5]"],tas=["sar/tas/*204[5-9]","sar/tas/*205[0-5]"])

months=['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
month_lengths=[31,28,31,30,31,30,31,31,30,31,30,31]
months_doy=[]
for i in range(len(month_lengths)):
    months_doy.extend([i]*month_lengths[i])

def load_geo(filename):
    """docstring for load_geo"""
    latnames=["lat","latitude","XLAT","XLAT_M"]
    lonnames=["lon","longitude","XLONG","XLONG_M"]
    lat=None
    for i in range(len(latnames)):
        try:
            lat=io.read_nc(filename,latnames[i]).data
            lon=io.read_nc(filename,lonnames[i]).data
        except:
            pass
    
    if (lat==None):
        raise KeyError
    
    if lon.max()>180:
        lon[lon>180]=lon[lon>180]-360
    
    if len(lat.shape)==1:
        lon,lat=np.meshgrid(lon,lat)
    if len(lat.shape)==3:
        lat=lat[0,...]
        lon=lon[0,...]
    return Bunch(lat=lat,lon=lon)
            

def load_monthly_ccsm(filename,geolut):
    """Load data from filename using geolut to interpolate spatially"""
    data=io.read_nc(filename,ccsm_var_to_read).data
    data=data[:12,:,:]*24*60*60*1000 # convert m/s to mm/day
    outputdata=np.zeros((12,geolut.shape[0],geolut.shape[1]))
    for i in range(4):
        y=geolut[:,:,i,0].astype('i')
        x=geolut[:,:,i,1].astype('i')
        outputdata+=np.float32(data[:,y,x]*geolut[np.newaxis,:,:,i,2])
    
    return outputdata

def doy2month(doy):
    """convert day of year array to month array"""
    return np.array([months_doy[int(thisdoy)] for thisdoy in doy])

def read_esg_data(filename,startyear,endyear):
    # read in daily data between start and end years and convert to monthly averaged data
    varname=esgvarname
    days=io.read_nc(filename,"time").data
    doy=days%365 #noleap calendar
    months=doy2month(doy)
    years=np.floor(days/365)
    ncdata=io.read_nc(filename,varname,returnNCvar=True)
    
    outputdata=np.zeros((12,ncdata.data.shape[1],ncdata.data.shape[2]))
    startpoint=np.where(years>=startyear)[0][0]
    endpoint=np.where(years<endyear)[0][-1]+1
    curdata=ncdata.data[startpoint:endpoint,:,:]
    ncdata.ncfile.close()

    for thismonth in range(12):
        curmonths=np.where(months[startpoint:endpoint]==thismonth)[0]
        current=curdata[curmonths,:,:]
        outputdata[thismonth,:,:]=current.mean(axis=0)
        
    if esgvarname=="pr":
        outputdata*=24*60*60 # convert mm/s to mm/day
    return outputdata

def load_monthly_esg_ccsm(filename,geolut):
    
    cur_data=read_esg_data(filename,1995,2006)
    fut_data=read_esg_data(filename,2045,2056)
    
    write_monthly_pgw(cur_data,fut_data,filename)
    # geographic interpolation
    current=np.zeros((12,geolut.shape[0],geolut.shape[1]))
    for i in range(4):
        y=geolut[:,:,i,0].astype('i')
        x=geolut[:,:,i,1].astype('i')
        current+=np.float32(cur_data[:,y,x]*geolut[np.newaxis,:,:,i,2])

    future=np.zeros((12,geolut.shape[0],geolut.shape[1]))
    for i in range(4):
        y=geolut[:,:,i,0].astype('i')
        x=geolut[:,:,i,1].astype('i')
        future+=np.float32(fut_data[:,y,x]*geolut[np.newaxis,:,:,i,2])
    
    
    return current,future

def load_ccsm(current,future):
    """Load CCSM mean monthly precipitation over the headwaters domain"""
    
    # ccsmgeo=load_geo(current)
    esggeo=load_geo(esgfile)
    wrfgeo=load_geo(geofile)
    # geoLUT=load_geoLUT(ccsmgeo.lat,ccsmgeo.lon,wrfgeo.lat,wrfgeo.lon)
    esggeoLUT=load_geoLUT(esggeo.lat,esggeo.lon,wrfgeo.lat,wrfgeo.lon)
    
    # this are the data used to compute the AR downscaled output
    # downloaded from ftp-esg.ucllnl.org/20c3m/atm/da/pr/ncar_ccsm3_0/run5/
    # concatenated with ftp-esg.ucllnl.org/sresa2/atm/da/pr/ncar_ccsm3_0/run5/
    cur_esg,fut_esg=load_monthly_esg_ccsm(esgfile,esggeoLUT)
    
    cur_data=cur_esg
    fut_data=fut_esg
    # these are the data from Kyoko that were used to generate the PGW simulations
    # they are *almost* identical (possibly off by a day due to the leap year / noleap calendar???)
    # cur_data=load_monthly_ccsm(current,geoLUT)
    # fut_data=load_monthly_ccsm(future,geoLUT)
    
    return Bunch(current=cur_data,future=fut_data), Bunch(current=cur_esg,future=fut_esg)

def load_monthly_wrf_precip(filesearch):
    """docstring for load_monthly_wrf"""
    files=glob.glob(filesearch)
    files.sort()
    files=files[:96]
    cur_data=io.read_nc(files[0],"RAINNC",returnNCvar=True)
    output_data=np.zeros((12,cur_data.data.shape[1],cur_data.data.shape[2]))
    last_data=cur_data.data[-1,...]*0
    cur_data.ncfile.close()
    
    for i in range(len(files)):
        curmonth=((i+9)%12)
        ncdata=io.read_nc(files[i],"RAINNC",returnNCvar=True)
        curdata=ncdata.data[-1,...]
        ncdata.ncfile.close()
        
        output_data[curmonth,...]+=(curdata-last_data)/month_lengths[curmonth]
        last_data=curdata
    
    output_data/=8.0
    return output_data

def load_monthly_wrf_tair(filesearch):
    """docstring for load_monthly_wrf"""
    files=glob.glob(filesearch)[0] #air T is all in one file
    data=io.read_nc(files,"T2").data[24:(-31*24)]

    output_data=np.zeros((12,data.shape[1],data.shape[2]))
    time0=datetime.datetime(2000,10,1,0,0,0)
    dt=datetime.timedelta(1/24.0)
    times=[time0+dt*i for i in range(data.shape[0])]

    newdelta=datetime.timedelta(0)
    # no leap years in this dataset, so bump up the timeoffset periodically
    for i in range(len(times)):
        if (times[i].month==2) and (times[i].day==29):
            newdelta+=datetime.timedelta(1)
        times[i]+=newdelta
        
    for i in range(data.shape[0]):
        curmonth=times[i].month-1
        output_data[curmonth,...]+=data[i,...]
    
    for i in range(12):
        output_data[i,...]/=month_lengths[i]*24
    
    return output_data

def load_monthly_wrf(filesearch):
    if esgvarname=="pr":
        return load_monthly_wrf_precip(filesearch[esgvarname])
    else:
        return load_monthly_wrf_tair(filesearch[esgvarname])

def load_wrf(current_search,future_search):
    current=load_monthly_wrf(current_search)
    future=load_monthly_wrf(future_search)
    return Bunch(current=current,future=future)
 
def load_monthly_ar(filesearch,geolut):
    """load monthly means for a given filesearch"""
    files=[]
    for search in filesearch:
        files.extend(glob.glob(search))
    files.sort()
    ncdata=io.read_nc(files[0],esgvarname,returnNCvar=True)
    outputdata=np.zeros((12,ncdata.data.shape[1],ncdata.data.shape[2]))
    ncdata.ncfile.close()
    
    for f in files:
        data=io.read_nc(f,esgvarname).data
        startpoint=0
        for i in range(12):
            endpoint=startpoint+month_lengths[i]
            outputdata[i,...]+=data[startpoint:endpoint,:,:].mean(axis=0)
            startpoint=endpoint
    outputdata/=len(files)
    
    data=np.zeros((12,geolut.shape[0],geolut.shape[1]))
    for i in range(4):
        y=geolut[:,:,i,0].astype('i')
        x=geolut[:,:,i,1].astype('i')
        data+=np.float32(outputdata[:,y,x]*geolut[np.newaxis,:,:,i,2])
    
    return data
            

def load_ar(current_search,future_search):
    """load monthly mean data for current and future datasets"""
    argeo=load_geo(glob.glob(current_search[0])[0])
    wrfgeo=load_geo(geofile)
    geoLUT=load_geoLUT(argeo.lat,argeo.lon,wrfgeo.lat,wrfgeo.lon)
    
    current=load_monthly_ar(current_search,geoLUT)
    future=load_monthly_ar(future_search,geoLUT)
    return Bunch(current=current,future=future)

def plot_maps(data,file_prefix):
    """Make an annual change map and a series of monthly change maps"""
    
    # set up for precip and temperature separately
    if esgvarname=="pr":
        curdata=(data.future-data.current).sum(axis=0)/12.0*365.25
        title="Change in Mean Annual Precipitation"
        vmin=-200; vmax=200
        mvmin=-2.5;mvmax=2.5
        cmap=plt.cm.seismic_r
    else:
        curdata=(data.future-data.current).mean(axis=0)
        title="Change in Mean Annual Temperature"
        vmin=-6; vmax=6
        mvmin=-5;mvmax=5
        cmap=plt.cm.seismic
        
    # annual map
    plt.figure();
    map_vis.vis(curdata,geo="subset",title=title,
                vmin=vmin,vmax=vmax,proj='lcc',
                cmap=cmap,colorbar=True,latstep=2.0,lonstep=4.0)
    plt.savefig(file_prefix+"_"+esgvarname+"_annual_map.png")
    plt.close()
    
    # Monthly plots
    plt.figure(figsize=(24,12))
    m=[]
    for i in range(12):
        plt.subplot(3,4,i+1)
        curdata=(data.future-data.current)[i,...]
        map_vis.vis(curdata,geo="subset",title=months[i],vmin=mvmin,vmax=mvmax,proj='lcc',
                    cmap=cmap,colorbar=True,latstep=2.0,lonstep=4.0,m=m)
        
        
    plt.savefig(file_prefix+"_"+esgvarname+"_monthly_maps.png")
    plt.close()
    
def plot_monthly_deltas(data,labels,prefix="",subset=None):
    """docstring for plot_monthly_deltas"""
    colors=["Blue","Red","Green","cyan","magenta"][:len(data)]
    
    if subset==None:
        for l,d,c in zip(labels,data,colors):
            plt.plot((d.future-d.current).mean(axis=1).mean(axis=1),label=l,color=c)
    else:
        for l,d in zip(labels,data):
            plt.plot((d.future-d.current)[:,subset[0]:subset[1],subset[2]:subset[3]].mean(axis=1).mean(axis=1),label=l)
        
    plt.legend()
    # plt.ylim(-1,1)
    plt.plot([0,11],[0,0],"--",color="k")
    plt.xlim(-0.5,11.5)
    plt.savefig(prefix+esgvarname+"_Monthly_deltas.png")
    plt.close()

def write_monthly_pgw(current,future,geofile):
    """write out a file containing the ratios future/current for each month and lat/lon data"""
    
    print("writing PGW results")
    if esgvarname=="pr":
        output_data=future/current
        outputfilename="PGW_ratio_file.nc"
        data_atts=Bunch(long_name="Ratio between future:current monthly precipitaton", units="mm/mm")
    else:
        output_data=future-current
        outputfilename="PGW_difference_file.nc"
        data_atts=Bunch(long_name="Difference between future-current monthly Temperature", units="K")
        
    lat=io.read_nc(geofile,"lat").data
    lon=io.read_nc(geofile,"lon").data
    time=np.arange(12)
    
    lat_atts =Bunch(long_name="latitude",  units="degrees")
    lon_atts =Bunch(long_name="longitude", units="degrees")
    time_atts=Bunch(long_name="Month of the year",units="month",description="0=January,11=December,etc.")
    
    datadims=("time","lat","lon")
    evars=[ Bunch(data=lat, name="lat", dtype="f",attributes=lat_atts,dims=("lat",)),
            Bunch(data=lon, name="lon", dtype="f",attributes=lon_atts,dims=("lon",)),
            Bunch(data=time,name="time",dtype="f",attributes=time_atts,dims=("time",)),
            ]
    
    io.write(outputfilename,output_data,dtype="f",varname="data",dims=datadims,attributes=data_atts,extravars=evars)
    

def main():
    """Make a series of figures comparing the climate change component from CCSM and WRF"""

    print("Loading CCSM")
    ccsm_data,esg_data=load_ccsm(ccsm_cur,ccsm_fut)
    print("Plotting CCSM")
    plot_maps(esg_data,"CCSM")

    print("Loading AR")
    ar_data=load_ar(ar_current_list[esgvarname],ar_future_list[esgvarname])
    print("Plotting AR")
    plot_maps(ar_data,"AR")
    
    print("Loading WRF")
    wrf_data=load_wrf(wrf_cur,wrf_fut)
    print("Plotting WRF")
    plot_maps(wrf_data,"WRF")
    
    
    print("Plotting Monthly comparison")
    labels=["WRF","CCSM","AR"]
    plot_monthly_deltas([wrf_data,esg_data,ar_data],labels=labels)
    plot_monthly_deltas([wrf_data,esg_data,ar_data],labels=labels,prefix="Mtns_",subset=[50,200,100,200])
    plot_monthly_deltas([wrf_data,esg_data,ar_data],labels=labels,prefix="Plains_",subset=[50,200,200,280])
    
    plot_monthly_deltas([wrf_data,esg_data],labels=labels,prefix="WRF_CCSM_")
    plot_monthly_deltas([wrf_data,esg_data],labels=labels,prefix="WRF_CCSM_Mtns_",subset=[50,200,100,200])
    plot_monthly_deltas([wrf_data,esg_data],labels=labels,prefix="WRF_CCSM_Plains_",subset=[50,200,200,280])

if __name__ == '__main__':
    main()

# from make_ccsm_comparison import *
# ccsm_data=load_ccsm(ccsm_cur,ccsm_fut)
# plot_maps(ccsm_data,"CCSM")
