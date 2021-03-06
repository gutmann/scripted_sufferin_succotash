#!/usr/bin/env python
import argparse
import traceback
import sys
import os
from glob import glob
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap

import mygis as myio

stats_location="/d2/gutmann/usbr/stat_data/DAILY/down/stats/"
onemm_loc="conusprecip_1mm/"
zeromm_loc="conusprecip_0mm/"
onemm_loc=""
zeromm_loc=onemm_loc
e0_loc="precip_0mm/"
e0_1mm_loc="precip_1mm/"

base_colormap=cm.jet_r
diff_colormap=cm.seismic_r

def load_monthly(filesearch,location=stats_location+onemm_loc,mask=None):
    files=glob(location+filesearch)
    files.sort()
    outputdata=np.zeros(12)
    for i,f in enumerate(files):
        data=np.ma.array(myio.read_nc(f).data,mask=mask)
        if data.max()>500:
            data[mask]=1e20
            data=np.ma.array(data,mask=(data>1000))
            # print(f)
        outputdata[i]=np.mean(data)
    return outputdata

def load_scales(filesearch,location=stats_location,mask=None):
    files=glob(location+filesearch)
    files.sort()
    files=np.array(files)[[0,3,2,1]]
    
    outputdata=np.zeros(4)
    for i,f in enumerate(files):
        data=myio.read_nc(f).data
        # extreme events have a dimension for 2yr,10yr,50yr,100yr... we want 50yr
        if len(data.shape)>2:
            data=data[2,...]
        
        if i==0:
            data=np.ma.array(data,mask=mask)
        else:
            data=np.ma.array(data,mask=data>2000)
        data=np.ma.array(data,mask=~np.isfinite(data))

        outputdata[i]=np.mean(data)
    return outputdata


def plot_v_season(data,colors,linestyles,markers,labels,yrange):
    
    xlabels=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    ax=plt.gca()
    
    for d,c,l,m,lab in zip(data,colors,linestyles,markers,labels):
        ax.plot(d,marker=m,color=c,linestyle=l,label=lab,linewidth=2,markersize=12)
    
    plt.ylim(yrange[0],yrange[1])
    
    ax.xaxis.set_ticks([])
    xaxis_location=yrange[0]-0.03*(yrange[1]-yrange[0])
    for i in range(len(xlabels)):
        ax.text(i-0.5,xaxis_location,xlabels[i],rotation=45)
        ax.plot([i,i],yrange,linestyle=":",color="black")
    plt.xlim(-0.5,12-0.5)
    
def plot_v_scale(data,colors,linestyles,markers,labels,yrange,label_x=True):
    
    xlabels=["Full Res","HUC 8","HUC 4", "HUC 2"]
    
    for d,c,l,m,lab in zip(data,colors,linestyles,markers,labels):
        plt.plot(d,marker=m,color=c,linestyle=l,label=lab,linewidth=2,markersize=12)
    
    plt.ylim(yrange[0],yrange[1])
    
    ax=plt.gca()
    ax.xaxis.set_ticks([])
    if label_x:
        xaxis_location=yrange[0]-0.03*(yrange[1]-yrange[0])
        for i in range(len(xlabels)):
            ax.text(i-0.4,xaxis_location,xlabels[i],rotation=45)
            
    plt.xlim(-0.5,4-0.5)


def map_vis(data,geo=[],title="",vmin=None,vmax=None,cmap=None,showcolorbar=True,
            latlabels=[1,0,0,0],lonlabels=[0,0,0,1],cbar_label=None):
    """Plot a map of data using the bounds in geo=[lllat,urlat,lllon,urlon]
    
    Optionally specify a map title, min and max value and colormap
    """
    # geo=[35,43,-113,-101]
    geo=[25.125,52.875,-124.75,-67]
    m = Basemap(projection='cyl',llcrnrlat=geo[0],urcrnrlat=geo[1],\
                llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
    
    mapimg=m.imshow(data,vmin=vmin,vmax=vmax,cmap=cmap)
    m.drawparallels(np.arange(25,55,5.),labels=latlabels,dashes=[1,4])
    m.drawmeridians(np.arange(-120,-65,10.),labels=lonlabels,dashes=[1,4])
    m.drawstates(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    
    if showcolorbar:
        cb=m.colorbar()
        if cbar_label!=None:
            cb.set_label(cbar_label)
        
    plt.title(title,x=0.05,ha="left")

def get_mask():
    d=myio.read_nc(stats_location+onemm_loc+"obs-maurer.125-pr_full_res_annual_MAP.nc").data
    mask=np.empty(d.shape,dtype=bool)
    mask[:]=False
    mask[d>1e10]=True
    return mask

#--------------------------------------------------------

def mean_annual_fig():
    
    files=["CA-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "obs-maurer.125-pr_full_res_annual_MAP.nc",
            "forcing-ncep-pr_full_res_annual_MAP.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ccsm-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "future/SAR-ccsm-pr-BC12km_full_res_annual_MAP.nc"
            ]
    pgw_start=5
    titles=["BCCA", "BCSDm", "Observations", "NCEP", "AR","AR-CCSM current","AR-PGW","AR-CCSM future"]
    letterlabels=[chr(letter)+") " for letter in range(ord('a'),ord('a')+len(titles))]
    cbar_label="Mean Annual Precip [mm/yr]"
    cbars=[False,True]*4
    mask=get_mask()
    
    plt.figure(figsize=(15,13),dpi=300)
    for i,f in enumerate(files):
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<(len(files)-2):lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]

        if i<pgw_start:
            data=np.ma.array(myio.read_nc(stats_location+onemm_loc+f).data,mask=mask)
        else:
            data=np.ma.array(myio.read_nc(f).data,mask=mask)
        plt.subplot(4,2,i+1)
        map_vis(data,title=letterlabels[i]+titles[i],cmap=base_colormap,vmin=0,vmax=2000,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
    plt.subplots_adjust(wspace = -0.1,hspace=0.15)
        
    plt.savefig("FIG2_Mean_Annual_precip.png",dpi=150)
    plt.close()

def bias_fig(month=None):
    
    files=["CA-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ccsm-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "future/SAR-ccsm-pr-BC12km_full_res_annual_MAP.nc"
            ]
    pgw_start=3
    titles=["BCCA - obs", "BCSDm - obs", "AR - obs","AR-CCSM current - obs","AR-PGW - obs","AR-CCSM future - obs"]
    letterlabels=[chr(letter)+") " for letter in range(ord('a'),ord('a')+len(titles))]
    cbar_label="Precipitation Bias [mm/yr]"
    cbars=[False,True]*3
    mask=get_mask()
    
    obsfile="obs-maurer.125-pr_full_res_annual_MAP.nc"
    biasrange=500
    if month!=None:
        obsfile=obsfile.replace("annual",month)
        biasrange=200
    obs=np.ma.array(myio.read_nc(stats_location+onemm_loc+obsfile).data,mask=mask)
    
    plt.figure(figsize=(15,10),dpi=300)
    for i,f in enumerate(files):
        if month!=None:
            f=f.replace("annual",month)
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<(len(files)-2):lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        if i>=pgw_start:
            data=np.ma.array(myio.read_nc(f).data,mask=mask)-obs
        else:
            data=np.ma.array(myio.read_nc(stats_location+onemm_loc+f).data,mask=mask)-obs
        plt.subplot(3,2,i+1)
        map_vis(data,title=letterlabels[i]+titles[i],cmap=diff_colormap,vmin=-biasrange,vmax=biasrange,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
    plt.subplots_adjust(wspace = -0.0001,hspace=0.15)
    if month==None:
        plt.savefig("FIG3_precip_bias.png",dpi=150)
    else:
        plt.savefig(month+"_FIG3_precip_bias.png",dpi=150)
    plt.close()

def agu_6panel_fig():
    """docstring for agu_6panel_fig"""
    files=["CA-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SD-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "../WRF-12km_annual_MAP.nc"]
    titles=["BCCA - obs", "BCSDd - obs", "AR - obs", "BCSDm - obs"," WRF - obs"]
    letterlabels=[chr(letter)+") " for letter in range(ord('a'),ord('a')+len(titles))]
    cbar_label="Precipitation Bias [mm/yr]"
    cbars=[False,True]*3
    # cbars.append(True)
    mask=get_mask()
    

    obsfile="obs-maurer.125-pr_full_res_annual_MAP.nc"
    obs=np.ma.array(myio.read_nc(stats_location+zeromm_loc+obsfile).data,mask=mask)
    
    plt.figure(figsize=(15,10),dpi=300)
    # plt.figure(figsize=(7,5),dpi=300)
    for i,f in enumerate(files):
        
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<3:lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        data=np.ma.array(myio.read_nc(stats_location+zeromm_loc+f).data,mask=mask)-obs
        plt.subplot(3,2,i+1)
        map_vis(data,title=letterlabels[i]+titles[i],cmap=diff_colormap,vmin=-500,vmax=500,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
        if i==4:
            plt.plot([-111,-111,-102.5,-102.5,-111],[35.5,42.5,42.5,35.5,35.5],color="green",linewidth=2)
    plt.subplots_adjust(wspace = -0.0001,hspace=0.15)
        
    plt.savefig("AGU_FIG_precip_bias.png",dpi=150)
    plt.close()

def changemap_fig():
    
    cbar_label="Precipitation Change [mm/yr]"
    mask=get_mask()

    o1="obs19c/obs-maurer.125-pr_full_res_annual_MAP.nc"
    o2="obs20c/obs-maurer.125-pr_full_res_annual_MAP.nc"
    n1="forcing19c/forcing-ncep-pr_full_res_annual_MAP.nc"
    n2="forcing20c/forcing-ncep-pr_full_res_annual_MAP.nc"
    
    obstitle="a) Observations (2000-2008) - (1979-1999)"
    nceptitle="b) NCEP (2000-2008) - (1979-1999)"
    
    plt.figure(figsize=(15,3.5),dpi=300)
        
    lonlabels=[0,0,0,1]
    latlabels=[1,0,0,0]
    
    
    obs1=np.ma.array(myio.read_nc(stats_location+o1).data,mask=mask)
    dobs=np.ma.array(myio.read_nc(stats_location+o2).data,mask=mask)-obs1
    plt.subplot(1,2,1)
    map_vis(dobs,title=obstitle,cmap=diff_colormap,vmin=-500,vmax=500,
            showcolorbar=True,cbar_label=cbar_label,
            latlabels=latlabels,lonlabels=lonlabels)
    
    latlabels=[0,0,0,0]
    ncep1=np.ma.array(myio.read_nc(stats_location+n1).data,mask=mask)
    dncep=np.ma.array(myio.read_nc(stats_location+n2).data,mask=mask)-ncep1
    plt.subplot(1,2,2)
    map_vis(dncep,title=nceptitle,cmap=diff_colormap,vmin=-500,vmax=500,
            showcolorbar=True,cbar_label=cbar_label,
            latlabels=latlabels,lonlabels=lonlabels)

    plt.subplots_adjust(wspace = -0.05)
        
    plt.savefig("FIG4_precip_change.png",dpi=150)
    plt.close()


def monthly_precip_fig():
    mask=get_mask()
    sar=load_monthly("SAR-ncep-pr-BC12km_full_res_month*_MAP.nc",mask=mask)
    sdd=load_monthly("SD-ncep-pr-BC12km_full_res_month*_MAP.nc",mask=mask)
    sdm=load_monthly("SDmon_c-ncep-pr-BC12km_full_res_month*_MAP.nc",mask=mask)
    ca=load_monthly("CA-ncep-pr-BC12km_full_res_month*_MAP.nc",mask=mask)
    obs=load_monthly("obs-maurer.125-pr_full_res_month*_MAP.nc",mask=mask)
    ncep=load_monthly("forcing-ncep-pr_full_res_month*_MAP.nc",mask=mask)

    pgw=load_monthly("SAR-ncep-pr-BC12km_full_res_month*_MAP.nc",location="./",mask=mask)
    ccsmc=load_monthly("SAR-ccsm-pr-BC12km_full_res_month*_MAP.nc",location="./",mask=mask)
    ccsmf=load_monthly("SAR-ccsm-pr-BC12km_full_res_month*_MAP.nc",location="future/",mask=mask)

    
    plt.figure(figsize=(7,5))
    
    plot_v_season([sar,sdd,sdm,ca,obs,ncep,pgw,ccsmc,ccsmf],
                colors=["g","darkred","red","b","k","k","g","g","g"],
                linestyles=["-","-","-","-","-","--",":","--","--"],
                markers=["None","None","None","None",".","None","None","None","^"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP","AR-PGW","AR-CCSMc","AR-CCSMf"],
                yrange=[30,150])
                
    plt.ylabel("mm")
    plt.title("Mean Monthly Precipitation")
    plt.legend(loc=1,ncol=3)
    plt.savefig("FIG5_monthly_precip.png",dpi=150)
    


def interannual_fig():
    
    files=["CA-ncep-pr-BC12km_full_res_annual_interannual.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_interannual.nc",
            "obs-maurer.125-pr_full_res_annual_interannual.nc",
            "forcing-ncep-pr_full_res_annual_interannual.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_interannual.nc",
            "SAR-ccsm-pr-BC12km_full_res_annual_interannual.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_interannual.nc",
            "future/SAR-ccsm-pr-BC12km_full_res_annual_interannual.nc"
            ]
    pgw_start=5
    titles=["BCCA", "BCSDm", "Observations", "NCEP", "AR","AR-CCSM current","AR-PGW","AR-CCSM future"]
    letterlabels=[chr(letter)+") " for letter in range(ord('a'),ord('a')+len(titles))]
    cbar_label="Interannual Variability [mm/yr]"
    cbars=[False,True]*4
    mask=get_mask()
    
    
    plt.figure(figsize=(15,13),dpi=300)
    for i,f in enumerate(files):
        
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<(len(files)-2):lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        if i>=pgw_start:
            data=np.ma.array(myio.read_nc(f).data,mask=mask)
        else:
            data=np.ma.array(myio.read_nc(stats_location+zeromm_loc+f).data,mask=mask)
        plt.subplot(4,2,i+1)
        map_vis(data,title=letterlabels[i]+titles[i],cmap=base_colormap,vmin=0,vmax=500,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
    plt.subplots_adjust(wspace = -0.1,hspace=0.15)
        
    plt.savefig("FIG6_Interannual_precip.png",dpi=150)
    plt.close()




def monthly_interannual_fig():
    mask=get_mask()
    sar=load_monthly("SAR-ncep-pr-BC12km_full_res_month*_interannual.nc",mask=mask)
    sdd=load_monthly("SD-ncep-pr-BC12km_full_res_month*_interannual.nc",mask=mask)
    sdm=load_monthly("SDmon_c-ncep-pr-BC12km_full_res_month*_interannual.nc",mask=mask)
    ca=load_monthly("CA-ncep-pr-BC12km_full_res_month*_interannual.nc",mask=mask)
    obs=load_monthly("obs-maurer.125-pr_full_res_month*_interannual.nc",mask=mask)
    ncep=load_monthly("forcing-ncep-pr_full_res_month*_interannual.nc",mask=mask)
    pgw=load_monthly("SAR-ncep-pr-BC12km_full_res_month*_interannual.nc",location="./",mask=mask)
    ccsmc=load_monthly("SAR-ccsm-pr-BC12km_full_res_month*_interannual.nc",location="./",mask=mask)
    ccsmf=load_monthly("SAR-ccsm-pr-BC12km_full_res_month*_interannual.nc",location="future/",mask=mask)
    
    plt.figure(figsize=(7,5))
    
    plot_v_season([sar,sdd,sdm,ca,obs,ncep,pgw,ccsmc,ccsmf],
                colors=["g","darkred","red","b","k","k","g","g","g"],
                linestyles=["-","-","-","-","-","--",":","--","--"],
                markers=["None","None","None","None",".","None","None","None","^"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP","AR-PGW","AR-CCSMc","AR-CCSMf"],
                yrange=[10,80])
                
    plt.ylabel("mm")
    plt.title("Interannual Variation (std.dev.)")
    plt.legend(loc=2,ncol=2)
    plt.savefig("FIG7_monthly_interannual.png",dpi=150)


def ext_scale_conus():
    mask=get_mask()
    sar =load_scales(stats_location+zeromm_loc+"SAR-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    sdd =load_scales(stats_location+zeromm_loc+"SD-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    sdm =load_scales(stats_location+zeromm_loc+"SDmon_c-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    ca  =load_scales(stats_location+zeromm_loc+"CA-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    obs =load_scales(stats_location+zeromm_loc+"obs-maurer.125-pr_*_annual_extremes_nday1.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,obs],
                colors=["g","darkred","red","b","k","k"],
                linestyles=["-","-","-","-","-","--"],
                markers=["None","None","None","None",".","None"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP"],
                yrange=[0,200])
                
    plt.ylabel("mm/day 50yr return")
    plt.title("CONUS")
    plt.legend(loc=1,ncol=2)
    
def ext_scale_sub():
    mask=None
    sar =load_scales(e0_loc+"SARe0-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    sdd =load_scales(e0_loc+"SDe0-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    sdm =load_scales(e0_loc+"SDmon-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    ca  =load_scales(e0_loc+"CA-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    cae  =load_scales(e0_loc+"CAe0-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    obs =load_scales(e0_loc+"obs-maurer.125-pr_*_annual_extremes_nday1.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,cae,obs],
                colors=["g","darkred","red","b","skyblue","k"],
                linestyles=["-","-","-","-","-","-"],
                markers=["None","None","None","None","None","."],
                labels=["AR","BCSDd","BCSDm","BCCA","BCCAr","Obs."],
                yrange=[0,120])
                
    # plt.ylabel("mm/day 50yr return")
    plt.title("Subdomain")
    plt.legend(loc=1,ncol=2)
    

def extreme_scaling_fig():
    plt.figure(figsize=(12,5))
    
    plt.subplot(1,2,1)
    ext_scale_conus()
    plt.subplot(1,2,2)
    ext_scale_sub()
    
    plt.savefig("FIG8_extreme_scaling.png",dpi=200)

def extremes_map_fig(month=None):
    
    files=["CA-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc",
            "obs-maurer.125-pr_full_res_annual_extremes_nday1.nc",
            "forcing-ncep-pr_full_res_annual_extremes_nday1.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc",
            "future/SAR-ccsm-pr-BC12km_full_res_annual_extremes_nday1.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc",
            "future/SAR-ccsm-pr-BC12km_full_res_annual_extremes_nday1.nc",
            ]
    pgw_start=5
    titles=["BCCA", "BCSDm", "Observations", "NCEP", "AR","AR-CCSM current","AR-PGW","AR-CCSM future"]
    letterlabels=[chr(letter)+") " for letter in range(ord('a'),ord('a')+len(titles))]
    cbar_label="50yr Extreme (mm)"
    cbars=[False,True]*4
    mask=get_mask()
    vmin=0
    vmax=500
    return_period=2 #0=2yr 1=10yr 2=50yr 3=100yr
    
    plt.figure(figsize=(15,13),dpi=300)
    for i,f in enumerate(files):
        if month!=None:
            f=f.replace("annual",month)
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<(len(files)-2):lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        if i<pgw_start:
            print(stats_location+onemm_loc+f)
            data=np.ma.array(myio.read_nc(stats_location+onemm_loc+f).data[return_period,...],mask=mask)
        else:
            print(f)
            data=np.ma.array(myio.read_nc(f).data[return_period,...],mask=mask)
        plt.subplot(4,2,i+1)
        map_vis(data,title=letterlabels[i]+titles[i],cmap=base_colormap,vmin=vmin,vmax=vmax,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
    plt.subplots_adjust(wspace = -0.1,hspace=0.15)
    
    if month==None:
        plt.savefig("FIG8.1_extremes_map.png",dpi=150)
    else:
        plt.savefig(month+"_FIG8.1_extremes_map.png",dpi=150)
    plt.close()

    

def wetfrac_map_fig(month=None):
    
    files=["CA-ncep-pr-BC12km_full_res_annual_wetfrac.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_wetfrac.nc",
            "obs-maurer.125-pr_full_res_annual_wetfrac.nc",
            "forcing-ncep-pr_full_res_annual_wetfrac.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_wetfrac.nc",
            "SAR-ccsm-pr-BC12km_full_res_annual_wetfrac.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_wetfrac.nc",
            "future/SAR-ccsm-pr-BC12km_full_res_annual_wetfrac.nc",
            ]
    pgw_start=5
    titles=["BCCA", "BCSDm", "Observations", "NCEP", "AR","AR-CCSM current","AR-PGW","AR-CCSM future"]
    letterlabels=[chr(letter)+") " for letter in range(ord('a'),ord('a')+len(titles))]
    cbar_label="Wet Day Fraction (0mm)"
    cbars=[False,True]*4
    mask=get_mask()

    
    plt.figure(figsize=(15,13),dpi=300)
    for i,f in enumerate(files):
        if month!=None:
            f=f.replace("annual",month)
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<(len(files)-2):lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        if i<pgw_start:
            data=np.ma.array(myio.read_nc(stats_location+onemm_loc+f).data,mask=mask)
        else:
            data=np.ma.array(myio.read_nc(f).data,mask=mask)
        plt.subplot(4,2,i+1)
        map_vis(data,title=letterlabels[i]+titles[i],cmap=base_colormap,vmin=0,vmax=1,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
    plt.subplots_adjust(wspace = -0.1,hspace=0.15)
    
    if month==None:
        plt.savefig("FIG9_Wetfrac_map.png",dpi=150)
    else:
        plt.savefig(month+"_FIG9_Wetfrac_map.png",dpi=150)
    plt.close()

def wet0_scale_conus():
    mask=get_mask()
    sar =load_scales(zeromm_loc+"SAR-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdd =load_scales(zeromm_loc+"SD-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdm =load_scales(zeromm_loc+"SDmon_c-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    ca  =load_scales(zeromm_loc+"CA-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    obs =load_scales(zeromm_loc+"obs-maurer.125-pr_*_annual_wetfrac.nc",mask=mask)
    try:
        ncep=load_scales("forcing-ncep-pr_*_annual_interannual.nc",mask=mask)
    except: 
        print("Missing ncep interannual as f(scale)")
        ncep=obs-10
    pgw =load_scales("SAR-ncep-pr-BC12km_*_annual_wetfrac.nc",location="./",mask=mask)
    ccsmc=load_scales("SAR-ccsm-pr-BC12km_*_annual_wetfrac.nc",location="./",mask=mask)
    ccsmf=load_scales("SAR-ccsm-pr-BC12km_*_annual_wetfrac.nc",location="future/",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,obs,ncep,pgw,ccsmc,ccsmf],
                colors=["g","darkred","red","b","k","k","g","g","g"],
                linestyles=["-","-","-","-","-","--",":","--","--"],
                markers=["None","None","None","None",".","None","None","None","^"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP","AR-PGW","AR-CCSMc","AR-CCSMf"],
                yrange=[0,1],label_x=False)
                
    plt.ylabel("Wet Day Fraction")
    plt.title("CONUS (0mm)")
    plt.legend(loc=4,ncol=2)
    
def wet0_scale_sub():
    mask=None
    sar =load_scales(e0_loc+"SARe0-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdd =load_scales(e0_loc+"SDe0-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdm =load_scales(e0_loc+"SDmon-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    ca  =load_scales(e0_loc+"CA-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    cae  =load_scales(e0_loc+"CAe0-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    obs =load_scales(e0_loc+"obs-maurer.125-pr_*_annual_wetfrac.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,cae,obs],
                colors=["g","darkred","red","b","skyblue","k"],
                linestyles=["-","-","-","-","-","-"],
                markers=["None","None","None","None","None","."],
                labels=["AR","BCSDd","BCSDm","BCCA","BCCAr","Obs."],
                yrange=[0,1],label_x=False)
                
    plt.title("Subdomain (0mm)")
    plt.legend(loc=4,ncol=2)
    
def wet1_scale_conus():
    mask=get_mask()
    sar =load_scales(onemm_loc+"SAR-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdd =load_scales(onemm_loc+"SD-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdm =load_scales(onemm_loc+"SDmon_c-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    ca  =load_scales(onemm_loc+"CA-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    obs =load_scales(onemm_loc+"obs-maurer.125-pr_*_annual_wetfrac.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,obs],
                colors=["g","darkred","red","b","k","k"],
                linestyles=["-","-","-","-","-","--"],
                markers=["None","None","None","None",".","None"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP"],
                yrange=[0.1,0.5])
                
    plt.ylabel("Wet Day Fraction")
    plt.title("CONUS (1mm)")
    # plt.legend(loc=4,ncol=2)
    
def wet1_scale_sub():
    mask=None
    sar =load_scales(e0_1mm_loc+"SARe0-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdd =load_scales(e0_1mm_loc+"SDe0-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdm =load_scales(e0_1mm_loc+"SDmon-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    ca  =load_scales(e0_1mm_loc+"CA-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    cae  =load_scales(e0_1mm_loc+"CAe0-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    obs =load_scales(e0_1mm_loc+"obs-maurer.125-pr_*_annual_wetfrac.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,cae,obs],
                colors=["g","darkred","red","b","skyblue","k"],
                linestyles=["-","-","-","-","-","-"],
                markers=["None","None","None","None","None","."],
                labels=["AR","BCSDd","BCSDm","BCCA","BCCAr","Obs."],
                yrange=[0.1,0.5])
                
    plt.title("Subdomain (1mm)")
    # plt.legend(loc=4,ncol=2)



def wetfrac_scale_fig():
    # plt.figure(figsize=(12,10))
    plt.figure(figsize=(7,6))
    
    # plt.subplot(2,2,1)
    wet0_scale_conus()
    # plt.subplot(2,2,2)
    # wet0_scale_sub()

    # plt.subplot(2,2,3)
    # wet1_scale_conus()
    # plt.subplot(2,2,4)
    # wet1_scale_sub()
    # 
    # plt.subplots_adjust(hspace=0.1)
    
    plt.savefig("FIG10_wetfrac_scaling.png",dpi=200)




def wetspell_dryspell_fig():
    mask=get_mask()
    sar=load_monthly("SAR-ncep-pr-BC12km_full_res_month*_wetspell.nc",mask=mask)
    sdd=load_monthly("SD-ncep-pr-BC12km_full_res_month*_wetspell.nc",mask=mask)
    sdm=load_monthly("SDmon_c-ncep-pr-BC12km_full_res_month*_wetspell.nc",mask=mask)
    ca=load_monthly("CA-ncep-pr-BC12km_full_res_month*_wetspell.nc",mask=mask)
    obs=load_monthly("obs-maurer.125-pr_full_res_month*_wetspell.nc",mask=mask)
    ncep=load_monthly("forcing-ncep-pr_full_res_month*_wetspell.nc",mask=mask)
    pgw =load_monthly("SAR-ncep-pr-BC12km_full_res_month*_wetspell.nc",location="./",mask=mask)
    ccsmc=load_monthly("SAR-ccsm-pr-BC12km_full_res_month*_wetspell.nc",location="./",mask=mask)
    ccsmf=load_monthly("SAR-ccsm-pr-BC12km_full_res_month*_wetspell.nc",location="future/",mask=mask)
    
    plt.figure(figsize=(14,5))
    
    plt.subplot(1,2,1)
    plot_v_season([sar,sdd,sdm,ca,obs,ncep,pgw,ccsmc,ccsmf],
                colors=["g","darkred","red","b","k","k","g","g","g"],
                linestyles=["-","-","-","-","-","--",":","--","--"],
                markers=["None","None","None","None",".","None","None","None","^"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP","AR-PGW","AR-CCSMc","AR-CCSMf"],
                yrange=[1.5,6])
                
    plt.ylabel("Days")
    plt.title("Wetspell Length")
    plt.legend(loc=2,ncol=1)


    sar=load_monthly("SAR-ncep-pr-BC12km_full_res_month*_dryspell.nc",mask=mask)
    sdd=load_monthly("SD-ncep-pr-BC12km_full_res_month*_dryspell.nc",mask=mask)
    sdm=load_monthly("SDmon_c-ncep-pr-BC12km_full_res_month*_dryspell.nc",mask=mask)
    ca=load_monthly("CA-ncep-pr-BC12km_full_res_month*_dryspell.nc",mask=mask)
    obs=load_monthly("obs-maurer.125-pr_full_res_month*_dryspell.nc",mask=mask)
    ncep=load_monthly("forcing-ncep-pr_full_res_month*_dryspell.nc",mask=mask)
    pgw =load_monthly("SAR-ncep-pr-BC12km_full_res_month*_dryspell.nc",location="./",mask=mask)
    ccsmc=load_monthly("SAR-ccsm-pr-BC12km_full_res_month*_dryspell.nc",location="./",mask=mask)
    ccsmf=load_monthly("SAR-ccsm-pr-BC12km_full_res_month*_dryspell.nc",location="./",mask=mask)
    
    plt.subplot(1,2,2)
    
    plot_v_season([sar,sdd,sdm,ca,obs,ncep,pgw,ccsmc,ccsmf],
                colors=["g","darkred","red","b","k","k","g","g","g"],
                linestyles=["-","-","-","-","-","--",":","--","--"],
                markers=["None","None","None","None",".","None","None","None","^"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP","AR-PGW","AR-CCSMc","AR-CCSMf"],
                yrange=[2,40])
                
    plt.ylabel("Days")
    plt.title("Dryspell Length")
    plt.legend(loc=2,ncol=1)
    
    plt.savefig("FIG11_monthly_spells.png",dpi=150)

def print_stats():
    from scipy.stats import ttest_ind as ttest
    
    mask=get_mask()
    sar=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SAR-ncep-pr-BC12km_full_res_annual_dryspell.nc").data,mask=mask)
    sdd=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SD-ncep-pr-BC12km_full_res_annual_dryspell.nc").data,mask=mask)
    sdm=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SDmon_c-ncep-pr-BC12km_full_res_annual_dryspell.nc").data,mask=mask)
    ca=np.ma.array(myio.read_nc(stats_location+onemm_loc+"CA-ncep-pr-BC12km_full_res_annual_dryspell.nc").data,mask=mask)
    obs=np.ma.array(myio.read_nc(stats_location+onemm_loc+"obs-maurer.125-pr_full_res_annual_dryspell.nc").data,mask=mask)
    print("AR    dryspell= "+str(sar.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sar[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDd dryspell= "+str(sdd.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdd[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDm dryspell= "+str(sdm.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdm[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCCA  dryspell= "+str( ca.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(ca[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("Obs   dryspell= "+str(obs.mean())[:5])
    print("")
    
    sar=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SAR-ncep-pr-BC12km_full_res_annual_wetspell.nc").data,mask=mask)
    sdd=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SD-ncep-pr-BC12km_full_res_annual_wetspell.nc").data,mask=mask)
    sdm=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SDmon_c-ncep-pr-BC12km_full_res_annual_wetspell.nc").data,mask=mask)
    ca=np.ma.array(myio.read_nc(stats_location+onemm_loc+"CA-ncep-pr-BC12km_full_res_annual_wetspell.nc").data,mask=mask)
    obs=np.ma.array(myio.read_nc(stats_location+onemm_loc+"obs-maurer.125-pr_full_res_annual_wetspell.nc").data,mask=mask)
    print("AR    wetspell= "+str(sar.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sar[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDd wetspell= "+str(sdd.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdd[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDm wetspell= "+str(sdm.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdm[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCCA  wetspell= "+str( ca.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(ca[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("Obs   wetspell= "+str(obs.mean())[:5])
    print("")
    
    sar=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SAR-ncep-pr-BC12km_full_res_annual_wetfrac.nc").data,mask=mask)
    sdd=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SD-ncep-pr-BC12km_full_res_annual_wetfrac.nc").data,mask=mask)
    sdm=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SDmon_c-ncep-pr-BC12km_full_res_annual_wetfrac.nc").data,mask=mask)
    ca=np.ma.array(myio.read_nc(stats_location+onemm_loc+"CA-ncep-pr-BC12km_full_res_annual_wetfrac.nc").data,mask=mask)
    obs=np.ma.array(myio.read_nc(stats_location+onemm_loc+"obs-maurer.125-pr_full_res_annual_wetfrac.nc").data,mask=mask)
    print("AR    wetfrac= "+str(sar.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sar[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDd wetfrac= "+str(sdd.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdd[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDm wetfrac= "+str(sdm.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdm[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCCA  wetfrac= "+str( ca.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(ca[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("Obs   wetfrac= "+str(obs.mean())[:5])
    print("")

    sar=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SAR-ncep-pr-BC12km_full_res_annual_MAP.nc").data,mask=mask)
    sdd=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SD-ncep-pr-BC12km_full_res_annual_MAP.nc").data,mask=mask)
    sdd.mask[~np.isfinite(sdd)]=True
    sdd.mask[sdd>1e5]=True
    sdm=np.ma.array(myio.read_nc(stats_location+onemm_loc+"SDmon_c-ncep-pr-BC12km_full_res_annual_MAP.nc").data,mask=mask)
    sdm.mask[~np.isfinite(sdm)]=True
    sdm.mask[sdm>1e5]=True
    ca=np.ma.array(myio.read_nc(stats_location+onemm_loc+"CA-ncep-pr-BC12km_full_res_annual_MAP.nc").data,mask=mask)
    obs=np.ma.array(myio.read_nc(stats_location+onemm_loc+"obs-maurer.125-pr_full_res_annual_MAP.nc").data,mask=mask)
    print("AR    MAP= "+str(sar.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sar[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDd MAP= "+str(sdd.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdd[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDm MAP= "+str(sdm.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdm[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCCA  MAP= "+str( ca.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(ca[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("Obs   MAP= "+str(obs.mean())[:5])
    print("")

    sar=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"SAR-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    sdd=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"SD-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    sdd.mask[~np.isfinite(sdd)]=True
    sdd.mask[sdd>1e5]=True
    sdm=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"SDmon_c-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    sdm.mask[~np.isfinite(sdm)]=True
    sdm.mask[sdm>1e5]=True
    ca=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"CA-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    obs=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"obs-maurer.125-pr_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    print("AR    extremes_nday1= "+str(sar.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sar[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDd extremes_nday1= "+str(sdd.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdd[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDm extremes_nday1= "+str(sdm.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(sdm[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCCA  extremes_nday1= "+str( ca.mean())[:5])
    print('  pvalue = {0}'.format(round(ttest(ca[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("Obs   extremes_nday1= "+str(obs.mean())[:5])
    print("")

    

def main():
    """Create figures for BAMS(?) synthesis downscaling paper"""

    # months=["month{:02}".format(i) for i in range(1,13)]
    # for month in months:
    #     print(month)
    #     # extremes_map_fig(month)
    #     bias_fig(month)
    #     wetfrac_map_fig(month)
    
    print("Extremes map")
    # extreme_scaling_fig()
    extremes_map_fig()

    #### hucmap_fig()
    # print("mean annual")
    mean_annual_fig()
    print("bias")
    bias_fig()
    # agu_6panel_fig()
    # changemap_fig()
    print("monthly")
    monthly_precip_fig()
    print("interannual")
    interannual_fig()
    print("monthly interannual")
    monthly_interannual_fig()
    print("wet frac map")
    wetfrac_map_fig()
    print("wet frac scale")
    # wetfrac_scale_fig()
    print("spell lengths")
    wetspell_dryspell_fig()
    #### geostats_fig()
    #### obs_wetfrac_fig()
    
    # print_stats()

if __name__ == '__main__':
    main()