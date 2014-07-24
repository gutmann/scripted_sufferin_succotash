#!/usr/bin/env python
import glob

import numpy as np
import matplotlib.pyplot as plt

from bunch import Bunch
import mygis
from stat_down import map_vis

wrfloc="/d5/gutmann/cc-downscaling-test/wrfoutput/T2m/"
sdloc="/d5/gutmann/cc-downscaling-test/SDdata/DAILY/"
geo_file=sdloc+"down/SAR/ccsm/pr/BCSAR_pr_12km_2000.nc"
ccsmgeo_file=sdloc+"model/ccsm/tas/tas.run5.ccsm.1999.nc"

sd_methods=["model/ccsm/tas/tas.run5.ccsm",
            "down/CA/ccsm/tas/BCCA_12km_CCSM_tas",
            "down/SAR/ccsm/tas/BCSAR_tas_12km",
            "down/SD/ccsm/tas/BCSD_12km_CCSM_tas",
            "down/SDmon/ccsm/tas/BCSDdisag_12km_CCSM_tave"]
# sd_methods=["model/ccsm/tas",
#             "down/CA/ccsm/tas/BCCA_12km_CCSM_tas"]
# clims=dict(MAP=(-300,300),wetfrac=[-0.1,0.1])
ndays_dict=dict(annual=365.25,
                month01=31,month02=28   ,month03=31,
                month04=30,month05=31   ,month06=30,
                month07=31,month08=31   ,month09=30,
                month10=31,month11=30   ,month12=31)

def get_time_bounds(timeperiod):
    """docstring for get_time_bounds"""
    if timeperiod=="annual":
        return (0,None)
        
    month=int(timeperiod[-2:])
    start=0
    for i in range(month-1):
        start+=ndays_dict["month{:02}".format(i+1)]
    
    stop=start+ndays_dict[timeperiod]
    if month==12:
        stop=None

    return (start,stop)

def load_time_period(files,varname,timeperiod):
    start,stop=get_time_bounds(timeperiod)
    print(start,stop,timeperiod)
    data=[]
    for f in files:
        curdata=mygis.read_nc(f,varname,returnNCvar=True)
        data.append(curdata.data[start:stop,:,:])
        curdata.ncfile.close()
    
    return np.concatenate(data)

def load_sd_data(stat,time,bounds=None):
    """load downscaling methods statistic for time"""
    output=[]
    means=[]
    sdvar="tas"
    for sd in sd_methods:
        #load current climate data from 1995-2005
        curfiles=glob.glob(sdloc+sd+"*199[5-9]*")
        curfiles.extend(glob.glob(sdloc+sd+"*200[0-5]*"))
        curfiles.sort()
        print(curfiles[0])
        current=load_time_period(curfiles,sdvar,time)
        
        #load future climate data from 2045-2055
        futfiles=glob.glob(sdloc+sd+"*204[5-9]*")
        futfiles.extend(glob.glob(sdloc+sd+"*205[0-5]*"))
        futfiles.sort()
        print(futfiles[0])
        future=load_time_period(futfiles,sdvar,time)
        
        #calculate the mean climate change over this time period
        data=future.mean(axis=0)-current.mean(axis=0)
        
        #subset the domain spatially if desired
        if bounds:
            if sd[:5]=="model":
                print("CCSM")
                means.append(data[bounds[0][0]:bounds[0][1],bounds[0][2]:bounds[0][3]].mean())
            else:
                print("Downscaled")
                means.append(data[bounds[1][0]:bounds[1][1],bounds[1][2]:bounds[1][3]].mean())
                print(means[-1])
        
        #save the name of the current method
        name=sd.split("/")[1]
        
        output.append(Bunch(data=data,name=name))
        
    return output,means



def load_wrf_data(stat,time):
    """Load WRF data and calculate the relevant statistic"""
    output_wrf_file=stat+"_"+time+".nc"
    try:
        current=mygis.read_nc(wrfloc+"current_"+output_wrf_file).data
        future=mygis.read_nc(wrfloc+"future_"+output_wrf_file).data
    except:
        print("    Reading raw data...")
        wrfvar="T2"
        if time=="annual":
            current=np.concatenate(mygis.read_files(wrfloc+"daily_NARR*",wrfvar))
            future=np.concatenate(mygis.read_files(wrfloc+"daily_PGW*",wrfvar))
        else:
            month=time[-2:]
            current=np.concatenate(mygis.read_files(wrfloc+"daily_NARR*"+month+".nc",wrfvar))
            future=np.concatenate(mygis.read_files(wrfloc+"daily_PGW*"+month+".nc",wrfvar))
        print("    Generating statistics")
        if stat=="tmean":
            current=current.mean(axis=0)
            future=future.mean(axis=0)
        
        print("    Writing stats for future reference")
        mygis.write(wrfloc+"current_"+output_wrf_file,current)
        mygis.write(wrfloc+"future_"+output_wrf_file,future)
        
    return future-current

    

def visualize(stat,time,sddata,wrfdata,geo,fig=None,m=None):
    """docstring for visualize"""
    
    # vmin,vmax=clims[stat]
    # clim=(vmin,vmax)
    clim=None
    if time=="annual":
        clim=(0,3)
    else:
        clim=(0,5)
    ccsm_geo=geo[0]
    sd_geo=geo[1]
    cmap=plt.cm.seismic
    cmap=plt.cm.jet
    
    if fig==None:
        fig=plt.figure(figsize=(10,10))
    else:
        fig.clf()
    plt.subplot(3,2,1)
    
    # first call to map_vis should be for WRF to set up the basemap
    if m==None:
        m=map_vis.vis(wrfdata,title="WRF",cmap=cmap,proj="lcc",clim=clim,
                        latstep=2.0,lonstep=5.0)
    else:
        map_vis.vis(wrfdata,title="WRF",cmap=cmap,proj="lcc",clim=clim,
                    m=m, latstep=2.0,lonstep=5.0)

    
    # now loop over remaining methods mapping each one reprojecting to the WRF grid
    for i,sd in enumerate(sddata):
        plt.subplot(3,2,i+2)
        if sd.name=="ccsm":
            map_vis.vis(sd.data,title=sd.name,cmap=cmap,clim=clim,proj="lcc",
                        m=m,reproject=True,lat=ccsm_geo.lat,lon=ccsm_geo.lon,
                        latstep=2.0,lonstep=5.0)
        else:
            map_vis.vis(sd.data,title=sd.name,cmap=cmap,clim=clim,
                        m=m,reproject=True,lat=sd_geo.lat,lon=sd_geo.lon,
                        latstep=2.0,lonstep=5.0)
    
    fig.savefig(time+"_"+stat+".png",dpi=200)
    #return the figure so it can be reused... would returning the basemap instance help too? 
    return fig,m

def get_bounds(geo):
    wrflatmin,wrflatmax,wrflonmin,wrflonmax=(34.19,43.65,-114.35,-99.64)
    miny=np.where(geo.lat[:,0]>wrflatmin)[0][0]
    maxy=np.where(geo.lat[:,0]>wrflatmax)[0][0]
    minx=np.where(geo.lon[0,:]>wrflonmin)[0][0]
    maxx=np.where(geo.lon[0,:]>wrflonmax)[0][0]
    
    return (miny,maxy,minx,maxx)
        

def plot_timeseries(timeseries,sdinfo,stat):
    plt.close()
    x=range(1,13)
    
    wrftime=[t[-1] for t in timeseries]
    for i in range(len(sdinfo)):
        sdtime=[t[i] for t in timeseries]
        if sdinfo[i].name=="ccsm":
            plt.plot(x,sdtime,label=sdinfo[i].name,linewidth=2)
        else:
            plt.plot(x,sdtime,label=sdinfo[i].name)
        
    plt.plot(x,wrftime,label="WRF",color="black",linewidth=2)
    plt.legend(ncol=2)
    
    plt.plot([0,13],[0,0],":",color="grey")
    
    plt.xlim(0.1,12.9)
    plt.xlabel("Month")
    plt.ylabel("Change in "+stat)
    
    plt.savefig(stat.replace(" ","_")+"_timeseries.png",dpi=100)
    

def main():
    """Plot a series of maps comparing different downscaling methods climate change portrayals"""
    
    patterns=["tmean"]
    times=["month{:02}".format(month+1) for month in range(12)]
    times.append("annual")
    # times=["annual"]
    
    timeseries=[]
    wftimeseries=[]
    
    geo=mygis.read_geo(geo_file)
    # geo.lat=geo.lat[:-1,:-1]
    # geo.lon=geo.lon[:-1,:-1]
    bounds=get_bounds(geo)
    
    cgeo=mygis.read_geo(ccsmgeo_file)
    cbounds=get_bounds(cgeo)
    print(cbounds)
    
    fig=None
    basemap=None
    timeonly=False
    for t in times:
        for p in patterns:
            print(t,p)
            print("Loading SD data")
            sd,means=load_sd_data(p,t,bounds=(cbounds,bounds))
            print("Loading WRF data")
            wrf=load_wrf_data(p,t)
            
            if t!="annual":
                means.append(wrf.mean())
                if p=="tmean":
                    timeseries.append(means)
            
            if not timeonly:
                print("visualizing...")
                fig,basemap=visualize(p,t,sd,wrf,(cgeo,geo),fig=fig) #m=basemap)

    plot_timeseries(timeseries,sd,"Mean Temperature")


if __name__ == '__main__':
    main()