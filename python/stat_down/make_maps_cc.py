#!/usr/bin/env python
import glob
import re

import numpy as np
import matplotlib.pyplot as plt

from bunch import Bunch
import mygis
from stat_down import map_vis

wrfloc="/d5/gutmann/cc-downscaling-test/wrfoutput/pr/"
sdloc="/d5/gutmann/cc-downscaling-test/SDdata/DAILY/down/"
geo_file=sdloc+"SAR/ccsm/pr/BCSAR_pr_12km_2000.nc"

sd_methods=["forcing-ccsm","CA-ccsm","SAR-ccsm","SD-ccsm","SDmon-ccsm"]
clims=dict(MAP=(-300,300),
           wetfrac=[-0.1,0.1],
           # dryspell = (-5,5),
           # wetspell = (-5,5),
           # mean_dtr = (-9999,9999),
           # mean_tave = (0,10),
           # mean_tmax = (0,10),
           # mean_tmin = (0,10)
           )
           # extremes_nday1 = (-9999,9999),
           # interannual = (-9999,9999),
           # frostdays = (-9999,9999),
           # growing_season = (-9999,9999),
           # interannual_dtr = (-9999,9999),
           # interannual_tave = (-9999,9999),
           # interannual_tmax = (-9999,9999),
           # interannual_tmin = (-9999,9999),
ndays_dict=dict(annual=365.25,
                month01=31,month02=28.25,month03=31,
                month04=30,month05=31   ,month06=30,
                month07=31,month08=31   ,month09=30,
                month10=31,month11=30   ,month12=31)

plot_labels=dict(CCSM="CCSM",
                      CA="BCCA",
                      SAR="AR",
                      SD="BCSDd",
                      SDmon="BCSDm")

line_plot_colors=dict(CCSM="grey",
                      CA="blue",
                      SAR="green",
                      SD="darkred",
                      SDmon="red")


def calc_wetfrac(data,threshold=0.0):
    """docstring for calc_wetfrac"""
    wetdays=np.zeros(data.shape)
    wetdays[data>threshold]=1
    return wetdays.mean(axis=0)

def calc_ndays(time):
    """docstring for calc_ndays"""
    return ndays_dict[time]
    
def load_sd_data(stat,time,bounds=None):
    """load downscaling methods statistic for time"""
    output=[]
    means=[]
    print("SD data loaded from pre-calculated summary files")
    for sd in sd_methods:
        try:
            curfile=glob.glob(sdloc+"current/"+sd+"*full_res_"+time+"_"+stat)[0]
        except IndexError as e:
            print(sdloc+"current/"+sd+"*full_res_"+time+"_"+stat)
            raise IndexError
        current=mygis.read_nc(curfile).data
        try:
            futfile=glob.glob(sdloc+"future/"+sd+"*full_res_"+time+"_"+stat)[0]
        except IndexError as e:
            print(sdloc+"future/"+sd+"*full_res_"+time+"_"+stat)
            raise e
        future=mygis.read_nc(futfile).data
        
        name=sd.split("-")[0]
        if name=="forcing":
            name="CCSM"
            if stat=="MAP":
                current*=86400
                future*=86400
        data=future-current
        
        if bounds:
            means.append(data[bounds[0]:bounds[1],bounds[2]:bounds[3]].mean())
        output.append(Bunch(data=data,name=name,label=plot_labels[name],color=line_plot_colors[name]))
        
    return output,means



def load_wrf_data(stat,time):
    """Load WRF data and calculate the relevant statistic"""
    output_wrf_file=stat+"_"+time+".nc"
    try:
        current=mygis.read_nc(wrfloc+"current_"+output_wrf_file).data
        future=mygis.read_nc(wrfloc+"future_"+output_wrf_file).data
        print("WRF data loaded from pre-calculated summary files")
    except:
        print("    Reading raw data...")
        if time=="annual":
            current=np.concatenate(mygis.read_files(wrfloc+"daily_NARR*"))
            future=np.concatenate(mygis.read_files(wrfloc+"daily_PGW*"))
        else:
            month=time[-2:]
            current=np.concatenate(mygis.read_files(wrfloc+"daily_NARR*"+month+".nc"))
            future=np.concatenate(mygis.read_files(wrfloc+"daily_PGW*"+month+".nc"))
        print("    Generating statistics")
        if stat=="MAP":
            ndays=calc_ndays(time)
            current=current.mean(axis=0)*ndays
            future=future.mean(axis=0)*ndays
        
        elif stat=="wetfrac":
            current=calc_wetfrac(current)
            future=calc_wetfrac(future)
        else:
            raise KeyError("stat not created for WRF:"+stat)
        print("    Writing stat for future reference")
        mygis.write(wrfloc+"current_"+output_wrf_file,current)
        mygis.write(wrfloc+"future_"+output_wrf_file,future)
        
    return future-current

def find_range(data):
    """docstring for find_range"""
    vmin=1e10
    vmax=-1e10
    for d in data:
        vmin=min(vmin,d.data.min())
        vmax=max(vmax,d.data.max())
    fullrange=max(abs(vmin),abs(vmax))
    
    return -fullrange,fullrange

def visualize(stat,time,sddata,wrfdata,geo,fig=None):
    """docstring for visualize"""
    
    vmin,vmax=clims[stat]
    if vmin==-9999:
        vmin,vmax=find_range(sddata)
    if time!="annual":
        if stat=="MAP":
            vmin/=5
            vmax/=5
        elif stat=="wetfrac":
            vmin*=2
            vmax*=2
    clim=(vmin,vmax)
    
    if fig==None:
        fig=plt.figure(figsize=(10,10))
    else:
        fig.clf()
    plt.subplot(3,2,1)
    if (stat=="MAP") or (stat=="wetfrac") or (stat=="wetspell"):
        cmap=plt.cm.seismic_r
    elif re.match("mean_t...",stat):
        cmap=plt.cm.jet
    else:
        cmap=plt.cm.seismic
        
        
    m=map_vis.vis(wrfdata,title="WRF",cmap=cmap,proj="lcc",clim=clim,
                    latstep=2.0,lonstep=5.0)
    
    for i,sd in enumerate(sddata):
        plt.subplot(3,2,i+2)
        map_vis.vis(sd.data,title=sd.label,cmap=cmap,clim=clim,
                    m=m,reproject=True,lat=geo.lat,lon=geo.lon,
                    latstep=2.0,lonstep=5.0)
    
    fig.savefig(time+"_"+stat+".png",dpi=300)
    # fig.savefig(time+"_"+stat+".pdf",dpi=300)
    return fig

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
    
    reorder=np.argsort([3,4,5,6,7,8,9,10,11,0,1,2])
    
    wrftime=np.array([t[-1] for t in timeseries])
    for i in range(len(sdinfo)):
        sdtime=np.array([t[i] for t in timeseries])
        if sdinfo[i].name=="CCSM":
            plt.plot(x,sdtime[reorder],label=sdinfo[i].label,linewidth=2,color=sdinfo[i].color)
        else:
            plt.plot(x,sdtime[reorder],label=sdinfo[i].label,color=sdinfo[i].color)
        
    plt.plot(x,wrftime[reorder],label="WRF",color="black",linewidth=2)
    plt.legend(ncol=2)
    
    plt.plot([0,13],[0,0],":",color="grey")
    
    plt.xlim(0.1,12.9)
    plt.xlabel("Month")
    plt.ylabel("Change in "+stat)
    xlabels=np.array(["J","F","M","A","M","J","J","A","S","O","N","D"])
    plt.xticks(x,xlabels[reorder])
    
    plt.savefig(stat.replace(" ","_")+"_timeseries.pdf",dpi=300)
    plt.savefig(stat.replace(" ","_")+"_timeseries.png",dpi=300)
    

def main():
    """Plot a series of maps comparing different downscaling methods climate change portrayals"""
    
    variable_names=["wetfrac","MAP"]
    variable_names=clims.keys()
    # variable_names=["MAP"]
    times=["month{:02}".format(month+1) for month in range(12)]
    times.append("annual")
    # times=["annual"]
    
    timeseries=[]
    wftimeseries=[]
    
    geo=mygis.read_geo(geo_file)
    geo.lat=geo.lat[:-1,:-1]
    geo.lon=geo.lon[:-1,:-1]
    bounds=get_bounds(geo)
    print(geo.lon.shape)
    print(bounds)
    
    fig=None
    timeonly=False
    for t in times:
        for v in variable_names:
            print(t,v)
            print("Loading SD data")
            sd,means=load_sd_data(v,t,bounds=bounds)
            print("Loading WRF data")
            try:
                wrf=load_wrf_data(v,t)
            except:
                wrf=load_wrf_data("MAP",t)
                
            
            if t!="annual":
                means.append(wrf.mean())
                if v=="MAP":
                    timeseries.append(means)
                elif v=="wetfrac":
                    wftimeseries.append(means)
            
            if not timeonly:
                print("visualizing...")
                fig=visualize(v,t,sd,wrf,geo,fig=fig)

    plot_timeseries(timeseries,sd,"Mean Precipitation")
    plot_timeseries(wftimeseries,sd,"Wet Day Fraction")


if __name__ == '__main__':
    main()