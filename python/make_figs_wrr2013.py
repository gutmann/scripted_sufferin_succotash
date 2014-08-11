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

# from stat_down import myio
import mygis as myio

stats_location="/d2/gutmann/usbr/stat_data/DAILY/down/stats-r1/"
onemm_loc="conusprecip_1mm/"
zeromm_loc="conusprecip_0mm/"
subdomain_names=["Southeast","Southwest","Northwest","Central Rockies"]
sub_locs=["stats-se/","stats-sw/","stats-pnw/","stats-cohw/"]
sub_locs_1mm=["stats-se_1mm/","stats-sw_1mm/","stats-pnw_1mm/","stats-cohw_1mm/"]
subdomain_names=subdomain_names[1:]
sub_locs=sub_locs[1:] #remove the SE domain
sub_locs_1mm=sub_locs_1mm[1:] #remove the SE domain
# e0_loc="precip_0mm/"
# e0_1mm_loc="precip_1mm/"

def load_monthly(filesearch,mask=None):
    files=glob(stats_location+onemm_loc+filesearch)
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

def load_scales(filesearch,mask=None):
    files=glob(stats_location+filesearch)
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

        if re.match(".*stats-pnw.*/.*huc2.*",f):
            # for the PNW subdomain, the other huc2s are too small a portion of the domain, so just use the PNW HUC2
            outputdata[i]=data[0,-2]
        else:
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

def get_mask(i=None):
    if i==None:
        d=myio.read_nc(stats_location+zeromm_loc+"obs-maurer.125-pr_full_res_annual_MAP.nc").data
        d+=myio.read_nc(stats_location+zeromm_loc+"SDmon*ncep-pr-BC12km_full_res_annual_MAP.nc").data
        catest=myio.read_nc(stats_location+zeromm_loc+"CA-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...]
    else:
        d=myio.read_nc(stats_location+sub_locs[i]+"obs-maurer.125-pr_full_res_annual_MAP.nc").data
        d+=myio.read_nc(stats_location+sub_locs[i]+"SDmon*ncep-pr-BC12km_full_res_annual_MAP.nc").data
        catest=myio.read_nc(stats_location+sub_locs[i]+"CA-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...]
    mask=np.empty(d.shape,dtype=bool)
    mask[:]=False
    mask[d>1e10]=True
    mask[~np.isfinite(catest)]=True
    return mask

#--------------------------------------------------------

def mean_annual_fig():
    
    files=["CA-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SD-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "obs-maurer.125-pr_full_res_annual_MAP.nc",
            "../forcing20c/forcing-ncep-pr_full_res_annual_MAP.nc"]
    titles=["a) BCCA", "b) BCSDd", "c) AR", "d) BCSDm", "e) Observations", "f) NCEP"]
    cbar_label="Mean Annual Precip [mm/yr]"
    cbars=[False,True]*3
    mask=get_mask()
    
    plt.figure(figsize=(15,10),dpi=300)
    for i,f in enumerate(files):
        
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<4:lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        data=np.ma.array(myio.read_nc(stats_location+zeromm_loc+f).data,mask=mask)
        plt.subplot(3,2,i+1)
        map_vis(data,title=titles[i],cmap=cm.jet,vmin=0,vmax=2000,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
    plt.subplots_adjust(wspace = -0.1,hspace=0.15)
        
    plt.savefig("FIG2_Mean_Annual_precip.png",dpi=150)
    plt.close()

def bias_fig():
    
    files=["CA-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SD-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_MAP.nc"]
    titles=["a) BCCA - obs", "b) BCSDd - obs", "c) AR - obs", "d) BCSDm - obs"]
    cbar_label="Precipitation Bias [mm/yr]"
    cbars=[False,True]*3
    mask=get_mask()
    
    obsfile="obs-maurer.125-pr_full_res_annual_MAP.nc"
    obs=np.ma.array(myio.read_nc(stats_location+zeromm_loc+obsfile).data,mask=mask)
    
    plt.figure(figsize=(15,7),dpi=300)
    for i,f in enumerate(files):
        
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<2:lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        data=np.ma.array(myio.read_nc(stats_location+zeromm_loc+f).data,mask=mask)-obs
        plt.subplot(2,2,i+1)
        map_vis(data,title=titles[i],cmap=cm.seismic,vmin=-500,vmax=500,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
    plt.subplots_adjust(wspace = -0.0001,hspace=0.15)
        
    plt.savefig("FIG3_precip_bias.png",dpi=150)
    plt.close()

def agu_6panel_fig():
    """docstring for agu_6panel_fig"""
    files=["CA-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SD-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_MAP.nc",
            "../WRF-12km_annual_MAP.nc"]
    titles=["a) BCCA - obs", "b) BCSDd - obs", "c) AR - obs", "d) BCSDm - obs"," e) WRF - obs"]
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
        map_vis(data,title=titles[i],cmap=cm.seismic,vmin=-500,vmax=500,
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
    map_vis(dobs,title=obstitle,cmap=cm.seismic,vmin=-500,vmax=500,
            showcolorbar=True,cbar_label=cbar_label,
            latlabels=latlabels,lonlabels=lonlabels)
    
    latlabels=[0,0,0,0]
    ncep1=np.ma.array(myio.read_nc(stats_location+n1).data,mask=mask)
    dncep=np.ma.array(myio.read_nc(stats_location+n2).data,mask=mask)-ncep1
    plt.subplot(1,2,2)
    map_vis(dncep,title=nceptitle,cmap=cm.seismic,vmin=-500,vmax=500,
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
    ncep=load_monthly("../forcing20c/forcing-ncep-pr_full_res_month*_MAP.nc",mask=mask)
    
    plt.figure(figsize=(7,5))
    
    plot_v_season([sar,sdd,sdm,ca,obs],
                colors=["g","darkred","red","b","k","k"],
                linestyles=["-","-","-","-","-","--"],
                markers=["None","None","None","None",".","None"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP"],
                yrange=[30,120])
                
    plt.ylabel("mm")
    plt.title("Mean Monthly Precipitation")
    plt.legend(loc=1,ncol=2)
    plt.savefig("FIG5_monthly_precip.png",dpi=150)
    


def interannual_fig():
    
    files=["CA-ncep-pr-BC12km_full_res_annual_interannual.nc",
            "SD-ncep-pr-BC12km_full_res_annual_interannual.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_interannual.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_interannual.nc",
            "obs-maurer.125-pr_full_res_annual_interannual.nc",
            "../forcing20c/forcing-ncep-pr_full_res_annual_interannual.nc"]
    titles=["a) BCCA", "b) BCSDd", "c) AR", "d) BCSDm", "e) Observations", "f) NCEP"]
    cbar_label="Interannual Variability [mm/yr]"
    cbars=[False,True]*3
    mask=get_mask()
    
    plt.figure(figsize=(15,10),dpi=300)
    for i,f in enumerate(files):
        
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<4:lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        data=np.ma.array(myio.read_nc(stats_location+zeromm_loc+f).data,mask=mask)
        plt.subplot(3,2,i+1)
        map_vis(data,title=titles[i],cmap=cm.jet,vmin=0,vmax=500,
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
    ncep=load_monthly("../forcing20c/forcing-ncep-pr_full_res_month*_interannual.nc",mask=mask)
    
    plt.figure(figsize=(7,5))
    
    plot_v_season([sar,sdd,sdm,ca,obs],
                colors=["g","darkred","red","b","k","k"],
                linestyles=["-","-","-","-","-","--"],
                markers=["None","None","None","None",".","None"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP"],
                yrange=[10,60])
                
    plt.ylabel("mm")
    plt.title("Interannual Variation (std.dev.)")
    plt.legend(loc=2,ncol=2)
    plt.savefig("FIG7_monthly_interannual.png",dpi=150)


def ext_scale_conus(xlabel=False):
    mask=get_mask()
    sar =load_scales(zeromm_loc+"SAR-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    sdd =load_scales(zeromm_loc+"SD-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    sdm =load_scales(zeromm_loc+"SDmon_c-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    ca  =load_scales(zeromm_loc+"CA-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    obs =load_scales(zeromm_loc+"obs-maurer.125-pr_*_annual_extremes_nday1.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,obs],
                colors=["g","darkred","red","b","k","k"],
                linestyles=["-","-","-","-","-","--"],
                markers=["None","None","None","None",".","None"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP"],
                yrange=[0,180],label_x=xlabel)
                
    plt.ylabel("mm/day 50yr return")
    plt.title("CONUS")
    # plt.legend(loc=1,ncol=2)
    
def ext_scale_sub(i,ylabel=False,xlabel=False,legend=False):
    mask=get_mask(i)
    sar =load_scales(sub_locs[i]+"SAR-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    sdd =load_scales(sub_locs[i]+"SD-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    sdm =load_scales(sub_locs[i]+"SDmon-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    ca  =load_scales(sub_locs[i]+"CA-ncep-pr-BC12km_*_annual_extremes_nday1.nc",mask=mask)
    cae =load_scales(sub_locs[i]+"BCCA_*-ncep-pr-12km_*_annual_extremes_nday1.nc",mask=mask)
    obs =load_scales(sub_locs[i]+"obs-maurer.125-pr_*_annual_extremes_nday1.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,cae,obs],
                colors=["g","darkred","red","b","skyblue","k"],
                linestyles=["-","-","-","-","-","-"],
                markers=["None","None","None","None","None","."],
                labels=["AR","BCSDd","BCSDm","BCCA","BCCAr","Obs."],
                yrange=[0,180],label_x=xlabel)
                
    if ylabel: plt.ylabel("mm/day 50yr return")
    plt.title(subdomain_names[i])
    if legend:plt.legend(loc=1,ncol=2)
    

def extreme_scaling_fig():
    plt.figure(figsize=(12,10))
    
    plt.subplot(2,2,1)
    ext_scale_conus()
    for i in range(3):
        plt.subplot(2,2,i+2)
        ext_scale_sub(i,ylabel=((i%2)==1),xlabel=i>0,legend=(i==1))
    
    plt.subplots_adjust(hspace=0.1)
    plt.savefig("FIG8_extreme_scaling.png",dpi=200)
    

def wetfrac_map_fig():
    
    files=["CA-ncep-pr-BC12km_full_res_annual_wetfrac.nc",
            "SD-ncep-pr-BC12km_full_res_annual_wetfrac.nc",
            "SAR-ncep-pr-BC12km_full_res_annual_wetfrac.nc",
            "SDmon_c-ncep-pr-BC12km_full_res_annual_wetfrac.nc",
            "obs-maurer.125-pr_full_res_annual_wetfrac.nc",
            "../forcing20c/forcing-ncep-pr_full_res_annual_wetfrac.nc"]
    titles=["a) BCCA", "b) BCSDd", "c) AR", "d) BCSDm", "e) Observations", "f) NCEP"]
    cbar_label="Wet Day Fraction (0mm)"
    cbars=[False,True]*3
    mask=get_mask()
    
    plt.figure(figsize=(15,10),dpi=300)
    for i,f in enumerate(files):
        
        lonlabels=[0,0,0,1]
        latlabels=[1,0,0,0]
        if i<4:lonlabels=[0,0,0,0]
        if (i%2)==1:latlabels=[0,0,0,0]
        
        data=np.ma.array(myio.read_nc(stats_location+zeromm_loc+f).data,mask=mask)
        plt.subplot(3,2,i+1)
        map_vis(data,title=titles[i],cmap=cm.jet,vmin=0,vmax=1,
                showcolorbar=cbars[i],cbar_label=cbar_label,
                latlabels=latlabels,lonlabels=lonlabels)
    plt.subplots_adjust(wspace = -0.1,hspace=0.15)
        
    plt.savefig("FIG9_Wetfrac_map.png",dpi=150)
    plt.close()

def wet0_scale_conus():
    mask=get_mask()
    sar =load_scales(zeromm_loc+"SAR-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdd =load_scales(zeromm_loc+"SD-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdm =load_scales(zeromm_loc+"SDmon_c-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    ca  =load_scales(zeromm_loc+"CA-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    obs =load_scales(zeromm_loc+"obs-maurer.125-pr_*_annual_wetfrac.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,obs],
                colors=["g","darkred","red","b","k","k"],
                linestyles=["-","-","-","-","-","--"],
                markers=["None","None","None","None",".","None"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP"],
                yrange=[0,1],label_x=False)
                
    plt.ylabel("Wet Day Fraction")
    plt.title("CONUS (0mm)")
    # plt.legend(loc=4,ncol=2)
    
def wet0_scale_sub(i,xlabel=False,legend=False):
    mask=get_mask(i)
    sar =load_scales(sub_locs[i]+"SAR-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdd =load_scales(sub_locs[i]+"SD-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdm =load_scales(sub_locs[i]+"SDmon-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    ca  =load_scales(sub_locs[i]+"CA-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    cae =load_scales(sub_locs[i]+"BCCA_*-ncep-pr-12km_*_annual_wetfrac.nc",mask=mask)
    obs =load_scales(sub_locs[i]+"obs-maurer.125-pr_*_annual_wetfrac.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,cae,obs],
                colors=["g","darkred","red","b","skyblue","k"],
                linestyles=["-","-","-","-","-","-"],
                markers=["None","None","None","None","None","."],
                labels=["AR","BCSDd","BCSDm","BCCA","BCCAr","Obs."],
                yrange=[0,1],label_x=xlabel)
                
    plt.ylabel("Wet Day Fraction")
    plt.title("{} (0mm)".format(subdomain_names[i]))
    if legend:plt.legend(loc=4,ncol=2)
    
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
                yrange=[0.1,0.5],label_x=False)
                
    # plt.ylabel("Wet Day Fraction")
    plt.title("CONUS (1mm)")
    # plt.legend(loc=4,ncol=2)
    
def wet1_scale_sub(i,xlabel=False):
    mask=get_mask(i)
    sar =load_scales(sub_locs_1mm[i]+"SAR-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdd =load_scales(sub_locs_1mm[i]+"SD-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    sdm =load_scales(sub_locs_1mm[i]+"SDmon_c-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    ca  =load_scales(sub_locs_1mm[i]+"CA-ncep-pr-BC12km_*_annual_wetfrac.nc",mask=mask)
    cae =load_scales(sub_locs_1mm[i]+"BCCA_*-ncep-pr-12km_*_annual_wetfrac.nc",mask=mask)
    obs =load_scales(sub_locs_1mm[i]+"obs-maurer.125-pr_*_annual_wetfrac.nc",mask=mask)
    
    plot_v_scale([sar,sdd,sdm,ca,cae,obs],
                colors=["g","darkred","red","b","skyblue","k"],
                linestyles=["-","-","-","-","-","-"],
                markers=["None","None","None","None","None","."],
                labels=["AR","BCSDd","BCSDm","BCCA","BCCAr","Obs."],
                yrange=[0.1,0.5],label_x=xlabel)
                
    plt.title("{} (1mm)".format(subdomain_names[i]))
    # plt.legend(loc=4,ncol=2)



def wetfrac_scale_fig():
    plt.figure(figsize=(12,20))
    
    nx=4
    ny=2
    plt.subplot(nx,ny,1)
    wet0_scale_conus()
    plt.subplot(nx,ny,3)
    wet0_scale_sub(0)
    plt.subplot(nx,ny,5)
    wet0_scale_sub(1)
    plt.subplot(nx,ny,7)
    wet0_scale_sub(2,xlabel=True,legend=True)

    plt.subplot(nx,ny,2)
    wet1_scale_conus()
    plt.subplot(nx,ny,4)
    wet1_scale_sub(0)
    plt.subplot(nx,ny,6)
    wet1_scale_sub(1)
    plt.subplot(nx,ny,8)
    wet1_scale_sub(2,xlabel=True)
    
    plt.subplots_adjust(hspace=0.1)
    
    plt.savefig("FIG10_wetfrac_scaling.png",dpi=200)




def wetspell_dryspell_fig():
    mask=get_mask()
    sar=load_monthly("SAR-ncep-pr-BC12km_full_res_month*_wetspell.nc",mask=mask)
    sdd=load_monthly("SD-ncep-pr-BC12km_full_res_month*_wetspell.nc",mask=mask)
    sdm=load_monthly("SDmon_c-ncep-pr-BC12km_full_res_month*_wetspell.nc",mask=mask)
    ca=load_monthly("CA-ncep-pr-BC12km_full_res_month*_wetspell.nc",mask=mask)
    obs=load_monthly("obs-maurer.125-pr_full_res_month*_wetspell.nc",mask=mask)
    
    plt.figure(figsize=(14,5))
    
    plt.subplot(1,2,1)
    plot_v_season([sar,sdd,sdm,ca,obs],
                colors=["g","darkred","red","b","k","k"],
                linestyles=["-","-","-","-","-","--"],
                markers=["None","None","None","None",".","None"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP"],
                yrange=[1,5])
                
    plt.ylabel("Days")
    plt.title("Wetspell Length")
    plt.legend(loc=2,ncol=1)


    sar=load_monthly("SAR-ncep-pr-BC12km_full_res_month*_dryspell.nc",mask=mask)
    sdd=load_monthly("SD-ncep-pr-BC12km_full_res_month*_dryspell.nc",mask=mask)
    sdm=load_monthly("SDmon_c-ncep-pr-BC12km_full_res_month*_dryspell.nc",mask=mask)
    ca=load_monthly("CA-ncep-pr-BC12km_full_res_month*_dryspell.nc",mask=mask)
    obs=load_monthly("obs-maurer.125-pr_full_res_month*_dryspell.nc",mask=mask)
    
    plt.subplot(1,2,2)
    
    plot_v_season([sar,sdd,sdm,ca,obs],
                colors=["g","darkred","red","b","k","k"],
                linestyles=["-","-","-","-","-","--"],
                markers=["None","None","None","None",".","None"],
                labels=["AR","BCSDd","BCSDm","BCCA","Obs.","NCEP"],
                yrange=[5,35])
                
    plt.ylabel("Days")
    plt.title("Dryspell Length")
    plt.legend(loc=2,ncol=1)
    
    plt.savefig("FIG11_monthly_spells.png",dpi=150)

def print_this_stat(statistic,i=None):
    from scipy.stats import ttest_ind as ttest
    if i==None:
        location=stats_location+zeromm_loc
        mask=get_mask()
    else:
        location=stats_location+sub_locs[i]
        mask=get_mask(i)
    subslice=slice(0,None,1)
    if statistic=="extremes_nday1":
        subslice=slice(2,3)
        mask=mask.reshape((1,mask.shape[0],mask.shape[1]))
    
    sar=np.ma.array(myio.read_nc(location+"SAR-ncep-pr-BC12km_full_res_annual_{}.nc".format(statistic)).data[subslice,...],mask=mask)
    sdd=np.ma.array(myio.read_nc(location+"SD-ncep-pr-BC12km_full_res_annual_{}.nc".format(statistic)).data[subslice,...],mask=mask)
    sdm=np.ma.array(myio.read_nc(location+"SDmon*-ncep-pr-BC12km_full_res_annual_{}.nc".format(statistic)).data[subslice,...],mask=mask)
    ca =np.ma.array(myio.read_nc(location+"CA-ncep-pr-BC12km_full_res_annual_{}.nc".format(statistic)).data[subslice,...],mask=mask)
    obs=np.ma.array(myio.read_nc(location+"obs-maurer.125-pr_full_res_annual_{}.nc".format(statistic)).data[subslice,...],mask=mask)
    
    
    obsmean=obs.mean()
    print("AR    {}= {}      error = {}".format(statistic,str(sar.mean())[:5],str(sar.mean()-obsmean)[:5]))
    print('  pvalue = {0}'.format(round(ttest(sar[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDd {}= {}      error = {}".format(statistic,str(sdd.mean())[:5],str(sdd.mean()-obsmean)[:5]))
    print('  pvalue = {0}'.format(round(ttest(sdd[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCSDm {}= {}      error = {}".format(statistic,str(sdm.mean())[:5],str(sdm.mean()-obsmean)[:5]))
    print('  pvalue = {0}'.format(round(ttest(sdm[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("BCCA  {}= {}      error = {}".format(statistic,str( ca.mean())[:5],str( ca.mean()-obsmean)[:5]))
    print('  pvalue = {0}'.format(round(ttest(ca[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    if i!=None:
        cae=np.ma.array(myio.read_nc(location+"BCCA_*-ncep-pr-12km_full_res_annual_{}.nc".format(statistic)).data[subslice,...],mask=mask)
        print("BCCAr {}= {}      error = {}".format(statistic,str(cae.mean())[:5],str(cae.mean()-obsmean)[:5]))
        print('  pvalue = {0}'.format(round(ttest(cae[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    print("Obs   {}= ".format(statistic)+str(obsmean)[:5])
    

def print_stats():
    from scipy.stats import ttest_ind as ttest
    
    mask=get_mask()
    # sar=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"SAR-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    # sdd=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"SD-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    # sdd.mask[~np.isfinite(sdd)]=True
    # sdd.mask[sdd>1e5]=True
    # sdm=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"SDmon_c-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    # sdm.mask[~np.isfinite(sdm)]=True
    # sdm.mask[sdm>1e5]=True
    # ca=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"CA-ncep-pr-BC12km_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    # ca.mask[~np.isfinite(ca)]=True
    # ca.mask[ca>1e5]=True
    # obs=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"obs-maurer.125-pr_full_res_annual_extremes_nday1.nc").data[2,...],mask=mask)
    # print("AR    extremes_nday1= "+str(sar.mean())[:5])
    # print('  pvalue = {0}'.format(round(ttest(sar[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    # print("BCSDd extremes_nday1= "+str(sdd.mean())[:5])
    # print('  pvalue = {0}'.format(round(ttest(sdd[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    # print("BCSDm extremes_nday1= "+str(sdm.mean())[:5])
    # print('  pvalue = {0}'.format(round(ttest(sdm[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    # print("BCCA  extremes_nday1= "+str( ca.mean())[:5])
    # print('  pvalue = {0}'.format(round(ttest(ca[mask==False], obs[mask==False])[1]*1000.0)/1000.0))
    # print("Obs   extremes_nday1= "+str(obs.mean())[:5])
    #
    # obs=np.ma.array(myio.read_nc(stats_location+zeromm_loc+"obs-uw.0625-pr_full_res_annual_extremes_nday1.nc").data[2,...])
    # mask=np.empty(obs.shape,dtype=bool)
    # mask[:]=False
    # mask[obs>1e10]=True
    # mask[~np.isfinite(obs)]=True
    # obs.mask=mask
    # print("Obs   extremes_nday1= "+str(obs.mean())[:5])

    print("")
    
    for stat in ["MAP","wetfrac","wetspell","dryspell","extremes_nday1"]:
        print("\n{}\n".format(stat))
        print_this_stat(stat)
        for i in range(len(sub_locs)):
            print("\n{}\n".format(subdomain_names[i]))
            print_this_stat(stat,i)
        
def huc_map():
    from stat_down import plot_hucs
    plot_hucs.main()

def main():
    """Create figures for Journal of Climate 2013 downscaling paper"""

    function_list=[
        # huc_map,
        mean_annual_fig, 
        bias_fig, 
        changemap_fig, 
        monthly_precip_fig, 
        interannual_fig, 
        monthly_interannual_fig, 
        extreme_scaling_fig, 
        wetfrac_map_fig, 
        wetfrac_scale_fig, 
        wetspell_dryspell_fig]#,
        # print_stats]
    
    for f in function_list:
        print(f.__name__)
        try:
            f()
        except Exception as e:
            print("ERROR: "+str(e))
    
    # huc_map()
    # mean_annual_fig()
    # bias_fig()
    # changemap_fig()
    # monthly_precip_fig()
    # interannual_fig()
    # monthly_interannual_fig()
    # extreme_scaling_fig()
    # wetfrac_map_fig()
    # wetfrac_scale_fig()
    # wetspell_dryspell_fig()
    #### geostats_fig()
    #### obs_wetfrac_fig()
    
    # print_stats()
        

if __name__ == '__main__':
    main()