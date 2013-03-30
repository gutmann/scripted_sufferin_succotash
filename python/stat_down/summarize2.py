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
from stat_down import myio as io
from mpl_toolkits.basemap import Basemap

parser= argparse.ArgumentParser(description='Summarize plots of various statistics. ')
parser.add_argument('-model',dest="models",nargs="?",action='store',
            default="ncep",help="model forcing to test ([ncep],narr)")
parser.add_argument('-res',dest="resolution",nargs="?",action='store',
            default="12km",help="resolution to run ([12km],6km)")
parser.add_argument('-exp',dest="experiment",nargs="?",action='store',
            default="e0",help="experiment to test ([e0],e1,conus)")
parser.add_argument('-pr',dest="pr_threshold",nargs="?",action='store',
            default="0",help="precipitation threshold used ([0],1)")
parser.add_argument ('--plotmaps', action='store_true',
        default=False, help='Make maps of all statistics and errors (takes a long time) [False]', dest='plotmaps')

args = parser.parse_args()

plotmaps=args.plotmaps
experiment=args.experiment
thisres=args.resolution
model_base=args.models
pr_threshold=int(args.pr_threshold)
# experiment="conus"
# experiment="e0"
# experiment="e1"

# thisres="6km"
# thisres="12km"

# model_base="ncep"
# model_base="narr"

if thisres=="12km":
    obs_base="obs-maurer.125-pr"
    if model_base=="ncep":
        SDbase="-ncep-pr-BC12km"
    else:
        SDbase="-narr-pr-BC12km"

if thisres=="6km":
    obs_base="obs-uw.0625-pr"
    if model_base=="ncep":
        SDbase="-ncep-pr-BC6km"
    else:
        SDbase="-narr-pr-BC6km"

ob2="../obs19c/"+obs_base

extremes=["extremes_nday1","extremes_nday2","extremes_nday3",
          "extremes_nday4","extremes_nday5"]
objectives=["dryspell","interannual","MAP","wetfrac","wetspell"]
objectives=["MAP","wetfrac","dryspell","interannual","wetspell"]
methods=["CAe0","CA","SARe0","SDe0","SDmon"]
scales=["full_res","huc8","huc4","huc2"]
times=["annual","month01","month02","month03","month04","month05",
       "month06","month07","month08","month09","month10","month11",
       "month12","season1","season2","season3","season4"]
if experiment=="conus":
    methods=["CA","SAR","SD","SDmon_c"]
    # scales=["full_res"]
    # times=["annual"]# ,"month01","month02","month03","month04","month05",
#            "month06","month07","month08","month09","month10","month11",
#            "month12","season1","season2","season3","season4"]
if experiment=="e1":
    methods=["CAe1","CA","SARe1","SDe1","SDmon"]
    # scales=["full_res"]
    # times=["annual"]# ,"month01","month02","month03","month04","month05",
#            "month06","month07","month08","month09","month10","month11",
#            "month12","season1","season2","season3","season4"]
        

all_mapscales=dict(dryspell=[0.,10.],interannual=[0.,400.],MAP=[0.,1400.],wetfrac=[0.0,1.0],wetspell=[0.,30.0])
if experiment=="conus":
    all_mapscales=dict(dryspell=[0.,10.],interannual=[0.,500.],MAP=[0.,2000.],wetfrac=[0.0,1.0],wetspell=[0.,30.0])

methodnames=dict(CAe0="BCCAr",CA="BCCA",SARe0="SAR",SDe0="BCSDd",SDmon="BCSDm",
                 CAe1="BCCAr",SARe1="SAR",SDe1="BCSDd",
                 SAR="SAR",SD="BCSDd",SDmon_c="BCSDm")

def map_vis(data,title=[],vmin=None,vmax=None,barlabel=None,cmap=None,outputdir=None):
    # print("visualizing : "+title)
    if cmap==None:
        cmap=cm.jet
    if len(data.shape)>2:
        plotdata=data[0,:,:]
    else:
        plotdata=data
    plt.clf()
    ny,nx=plotdata.shape
    geo=[35,43,-113,-101]
    if experiment=="conus":geo=[25,52.7,-124.7,-67]
    m = Basemap(projection='cyl',llcrnrlat=geo[0],urcrnrlat=geo[1],\
                llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
    
    mapimg=m.imshow(plotdata,vmin=vmin,vmax=vmax,cmap=cmap)
    if experiment=="conus":
        m.drawparallels(np.arange(25,55,5.),labels=[1,0,0,0],dashes=[1,4])
        m.drawmeridians(np.arange(-120,-65,10.),labels=[0,0,0,1],dashes=[1,4])
        m.drawstates(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawcoastlines(linewidth=0.5)
    else:
        m.drawparallels(np.arange(36,43,2.),labels=[1,0,0,0],dashes=[1,4])
        m.drawmeridians(np.arange(-112,-103,4.),labels=[0,0,0,1],dashes=[1,4])
        m.drawstates(linewidth=1.5)
    cbar=m.colorbar()
    if barlabel:
        cbar.set_label(barlabel)
    plt.title(" ".join(title))
    if outputdir:
        plt.savefig(outputdir+"_".join(title)+'_map.png')
    else:
        plt.savefig("_".join(title)+'_map.png')

def titleize(name):
    namelist=list(name)
    namelist[0]=namelist[0].upper()
    return "".join(namelist)

def plot_vscale(data,plotname,yrange=None,plotcol=0,legend_loc=2,legend_ncol=1,obs=None,obs2=None,ylabel="Bias"):
    if experiment=="conus":
        sdm=3;sdd=2;sar=1;cac=0
        methods=["CA","SAR","SD","SDmon_c"]
        varpos=[sdm,sdd,sar,cac]
        names=["BCSDm","BCSDd","SAR","BCCA"]
        colors=["red","darkred","green","blue"]
    else:
        sdm=4;sdd=3;sar=2;cac=1;car=0
        varpos=[sdm,sdd,sar,cac,car]
        names=["BCSDm","BCSDd","SAR","BCCA","BCCAr"]
        colors=["red","darkred","green","blue","skyblue"]
    nscales=len(scales)
    
    plt.clf();
    ax=plt.gca()
    for v,n,c in zip(varpos,names,colors):
        ax.plot(data[v*nscales:(v+1)*nscales,plotcol],label=n,color=c,linewidth=2)
    if obs!=None:
        ax.plot(obs[0:nscales,plotcol],'-o',label="obs",color="k",linewidth=3)
    if obs2!=None:
        ax.plot(obs2[0:nscales,plotcol],'-x',label="obs2",color="k",linewidth=3)
    
    ax.xaxis.set_ticks([])
    if yrange!=None:
        xaxis_location=yrange[0]-0.03*(yrange[1]-yrange[0])
    else:
        xaxis_location=data.min()
    for i in range(nscales):
        ax.text(i-0.1,xaxis_location,scales[i].replace("huc","huc ").replace("_"," "),rotation=45)
    if yrange!=None:
        plt.ylim(yrange[0],yrange[1])
    plt.xlim(-0.5,nscales-0.5)
    plt.ylabel(ylabel)
    ax.plot([-0.5,nscales-0.5],[0,0],linestyle=":",color="black")
    plt.title(titleize(plotname))
    plt.legend(loc=legend_loc,ncol=legend_ncol)
    plt.draw()
    plotname="_".join([plotname,experiment,model_base,thisres])
    plt.savefig(plotname.replace(" ","_")+"_comp_scale.png")
    
def plot_vseason(data,plotname,yrange=None,plotrow=0,multiplier=1.0,legend_ncol=1,obs=None,obs2=None,ylabel="Bias"):
    if experiment=="conus":
        sdm=3;sdd=2;sar=1;cac=0
        methods=["CA","SAR","SD","SDmon_c"]
        varpos=[sdm,sdd,sar,cac]
        names=["BCSDm","BCSDd","SAR","BCCA"]
        colors=["red","darkred","green","blue"]
    else:
        sdm=4;sdd=3;sar=2;cac=1;car=0
        varpos=[sdm,sdd,sar,cac,car]
        names=["BCSDm","BCSDd","SAR","BCCA","BCCAr"]
        # colors=["red","magenta","green","blue","cyan"]
        colors=["red","darkred","green","blue","skyblue"]
    nscales=len(scales)
    xlabels=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    
    plt.clf();
    ax=plt.gca()
    for v,n,c in zip(varpos,names,colors):
        annual_value=str(data[v*nscales+plotrow,0]).split(".")
        if len(annual_value)>1:
            annual_value[1]=annual_value[1][:2]
        annual_value=".".join(annual_value)
        ax.plot(data[v*nscales+plotrow,1:13]/multiplier,label=n+" yr="+annual_value,color=c,linewidth=2)
    if obs!=None:
        annual_value=str(obs[0+plotrow,0]).split(".")
        if len(annual_value)>1:
            annual_value[1]=annual_value[1][:2]
        annual_value=".".join(annual_value)
        ax.plot(obs[0+plotrow,1:13]/multiplier,'-o',label="obs yr="+annual_value,color="k",linewidth=3)
    if obs2!=None:
        annual_value=str(obs2[0+plotrow,0]).split(".")
        if len(annual_value)>1:
            annual_value[1]=annual_value[1][:2]
        annual_value=".".join(annual_value)
        ax.plot(obs2[0+plotrow,1:13]/multiplier,'-x',label="obs2 yr="+annual_value,color="k",linewidth=3)
    
    ax.xaxis.set_ticks([])
    if yrange!=None:
        xaxis_location=yrange[0]-0.03*(yrange[1]-yrange[0])
    else:
        xaxis_location=data.min()
    for i in range(len(xlabels)):
        ax.text(i-0.1,xaxis_location,xlabels[i],rotation=45)
    if yrange!=None:
        plt.ylim(yrange[0],yrange[1])
    plt.xlim(-0.5,12-0.5)
    plt.ylabel(ylabel)
    ax.plot([-0.5,12-0.5],[0,0],linestyle=":",color="black")
    plt.title(titleize(plotname))
    plt.legend(loc=2,ncol=legend_ncol)
    plt.draw()
    plotname="_".join([plotname,experiment,model_base,thisres])
    plt.savefig(plotname.replace(" ","_")+"_comp_season.png")


def make_plot(data,minmax,name,obj,cmap=None):
    ny,nx=data.shape
    subfs=8
    xlabels=["Ann","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","DJF","MAM"," JJA","SON"]
    
    plt.clf()
    plt.imshow(data,vmin=minmax[0],vmax=minmax[1],cmap=cmap)
    
    # annotations and decorations
    for i in np.arange(-0.5,ny,len(scales)):plt.plot([-3.5,nx-0.5],[i,i],color='k',linewidth=2)
    plt.plot([0.5,0.5],[-1.5,ny-0.5],color='k',linewidth=2)
    plt.plot([12.5,12.5],[-1.5,ny-0.5],color='k',linewidth=2)
    for i in range(len(methods)):
        plt.text(-5,4*i+1.25,methodnames[methods[i]])
    for i in range(ny):
        plt.text(-2.2,i-0.25,scales[i%len(scales)].replace("_"," ").replace("huc","huc "),fontsize=subfs)
    for i in range(nx):
        plt.text(i-0.5,-1,xlabels[i],fontsize=subfs,rotation=45)
        
    ax=plt.gca()
    for i in range((len(methods))):
        ypos=4*i-0.5
        ax.annotate("",xy=(0, ypos), xycoords='data',
                    xytext=(-3.5,ypos), textcoords='data',
                    arrowprops=dict(arrowstyle="-",color='k'))
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    plt.xlim(-0.5,nx-0.5)
    plt.ylim(-0.5,ny-0.5)
    plt.title(titleize(name)+" Error "+titleize(obj))
    plt.colorbar()
    plt.draw()
    name="_".join([name,experiment,model_base,thisres])
    plt.savefig(obj+"_"+name+".png")


def extreme_plot(data,minmax,name,cmap=None):
    plt.clf()
    ny,nx=data.shape
    subfs=8
    xyearlabels=["  2yr"," 10yr"," 50yr","100yr"]
    
    plt.imshow(data,vmin=minmax[0],vmax=minmax[1],cmap=cmap)
    
    #annotations and decorations
    for i in np.arange(-0.5,ny,len(scales)):plt.plot([-3.5,nx-0.5],[i,i],color='k',linewidth=2)
    for i in range(len(methods)):
        plt.text(-5,4*i+1.25,methodnames[methods[i]])
    for i in range(ny):
        plt.text(-2.2,i-0.25,scales[i%len(scales)].replace("_"," ").replace("huc","huc "),fontsize=subfs)
    for i in range(len(extremes)):
        plt.plot([4*i-0.5,4*i-0.5],[-2,ny-0.5],color="k",linewidth=2)
        plt.text(4*i+0.25,-2.5,"ndays="+extremes[i][-1:],fontsize=subfs+1)
    for i in range(nx):
        plt.text(i-0.45,-1,xyearlabels[i%len(xyearlabels)],fontsize=subfs,rotation=45)
    ax=plt.gca()
    for i in range((len(methods))):
        ypos=4*i-0.5
        ax.annotate("",xy=(0, ypos), xycoords='data',
                    xytext=(-3.5,ypos), textcoords='data',
                    arrowprops=dict(arrowstyle="-",color='k'))
    for i in range((len(extremes))):
        xpos=4*i-0.5
        ax.annotate("",xy=(xpos, 0), xycoords='data',
                    xytext=(xpos,-3.5), textcoords='data',
                    arrowprops=dict(arrowstyle="-"),color='k')
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    plt.xlim(-0.5,nx-0.5)
    plt.ylim(-0.5,ny-0.5)
    plt.title("Fractional "+titleize(name)+" Error in Extremes")
    plt.colorbar()
    plt.draw()
    name="_".join([name,experiment,model_base,thisres])
    plt.savefig("extremes_"+name+".png")
    

def plot_histograms(obbase,sdbase,scale,time):
    obj="histogram"
    methods=["CAe0","CA","SARe0","SDe0","SDmon"]
    names=dict(SDmon="BCSDm",SDe0="BCSDd",SARe0="SAR",CA="BCCA",CAe0="BCCAr")
    colors=dict(SDmon="red",SDe0="darkred",SARe0="green",CA="blue",CAe0="skyblue")
    if experiment=="e1":
        methods=["CAe1","CA","SARe1","SDe1","SDmon"]
        names=dict(SDmon="BCSDm",SDe1="BCSDd",SARe1="SAR",CA="BCCA",CAe1="BCCAr")
        colors=dict(SDmon="red",SDe1="darkred",SARe1="green",CA="blue",CAe1="skyblue")
    if experiment=="conus":
        methods=["CA","SAR","SD","SDmon_c"]
        names=dict(SDmon="BCSDm",SD="BCSDd",SAR="SAR",CA="BCCA",CAe1="BCCAr",SDmon_c="BCSDm")
        colors=dict(SDmon="red",SD="darkred",SAR="green",CA="blue",CAe1="skyblue",SDmon_c="red")

    obs=io.read_nc("_".join([obbase,scale,time,obj])+".nc").data
    plt.clf()
    sd=[]
    for i in range(len(methods)):
        sd.append(io.read_nc("_".join([methods[i]+sdbase,scale,time,obj])+".nc").data)
        plt.plot(sd[-1][1,:],sd[-1][0,:],color=colors[methods[i]],label=names[methods[i]])
        
    plt.plot(obs[1,:],obs[0,:],'o-',linewidth=2,color="k",label="obs")
    plt.ylabel("n")
    plt.xlabel("Daily Precipitation [mm]")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylim(1,1E6)
    if experiment=="conus":
        plt.ylim(1,1E7)
    if re.match(".*huc2.*",scale):
        plt.ylim(1,3e3)
    if re.match(".*huc4.*",scale):
        plt.ylim(1,1e4)
    if re.match(".*huc8.*",scale):
        plt.ylim(1,1e5)
    plt.xlim(1E-2,400)
    plt.legend(loc=3,ncol=2)
    name="_".join(["hists",scale,time,experiment,model_base,thisres])
    plt.savefig(name+".png")
    
    
    
if __name__ == '__main__':
    
    example_data=io.read_nc(glob("_".join([obs_base,"full_res","annual","MAP"])+".nc")[0]).data
    maskit=((example_data>8000)|(example_data<0)|(~np.isfinite(example_data)))
    nx=len(times)
    ny=len(methods)*len(scales)
    if experiment=="conus":
        ny=5*len(scales)
        nx=12+4+1
    # calculating maps for each statistics (objective function)
    for obj in objectives:
        fullobs=np.zeros((ny,nx))
        fullobs2=np.zeros((ny,nx))
        RMS=np.zeros((ny,nx))
        bias=np.zeros((ny,nx))
        for j in range(len(scales)):
            for t in range(len(times)):
                mapscale=all_mapscales[obj]
                if (t>0) and ((obj=="MAP") or (obj=="interannual")):
                    mapscale=np.array(mapscale)/5
                obs=io.read_nc("_".join([obs_base,scales[j],times[t],obj])+".nc").data
                try:
                    obs2=io.read_nc("_".join([ob2,scales[j],times[t],obj])+".nc").data
                except Exception as e:
                    obs2=obs
                    print(e)
                if scales[j]=="full_res":
                    obs=np.ma.array(obs,mask=maskit)
                    obs2=np.ma.array(obs2,mask=maskit)
                else:
                    if ((obj=="MAP") or (obj=="interannual")):
                        obs=np.ma.array(obs,mask=obs>5000)
                        obs2=np.ma.array(obs2,mask=obs2>5000)
                    else:
                        obs=np.ma.array(obs,mask=obs>2900)
                        obs2=np.ma.array(obs2,mask=obs2>2900)
                        # print(obs_base,scales[j],times[t],obj,obs.max())
                        # obs=np.ma.array(obs,mask=obs>50)
                        
                if plotmaps:
                    map_vis(obs,[obs_base,scales[j],times[t],obj],vmin=mapscale[0],vmax=mapscale[1])
                if (obj=="MAP") and (times[t]=="annual"):
                    plot_histograms(obs_base,SDbase,scales[j],times[t])
                for i in range(len(methods)):
                    try:
                        sd=io.read_nc("_".join([methods[i]+SDbase,scales[j],times[t],obj])+".nc").data
                        for curdim in range(len(sd.shape)):
                            if sd.shape[curdim]!=obs.shape[curdim]:
                                if curdim==0:
                                    if obs.shape[curdim]<sd.shape[curdim]:
                                        sd=sd[:obs.shape[curdim],:]
                                    else:
                                        obs=obs[:sd.shape[curdim],:]
                                        obs2=obs2[:sd.shape[curdim],:]
                                if curdim==1:
                                    if obs.shape[curdim]<sd.shape[curdim]:
                                        sd=sd[:,:obs.shape[curdim]]
                                    else:
                                        obs=obs[:,:sd.shape[curdim]]
                                        obs2=obs2[:,:sd.shape[curdim]]
                        if scales[j]=="full_res":
                            sd=np.ma.array(sd,mask=((~np.isfinite(sd))|maskit|(sd>6000)))
                        else:
                            if ((obj=="MAP") or (obj=="interannual")):
                                sd=np.ma.array(sd,mask=((~np.isfinite(sd))|(obs>5000)|(sd>5000)))
                            else:
                                sd=np.ma.array(sd,mask=((~np.isfinite(sd))|(obs>2900)|(sd>5000)))
                                # print(methods[i],SDbase,scales[j],times[t],obj,sd.max())
                                # sd=np.ma.array(sd,mask=sd>50)
                                # print(methods[i],SDbase,scales[j],times[t],obj,sd.max())
                        if plotmaps:
                            map_vis(sd,[methods[i]+SDbase,scales[j],times[t],obj],vmin=mapscale[0],vmax=mapscale[1])
                        err=sd-obs
                        if plotmaps:
                            map_vis(err,['error',methods[i]+SDbase,scales[j],times[t],obj],
                                    cmap=cm.seismic,vmin=-mapscale[1]/2.0,vmax=mapscale[1]/2.0)
                        good=np.where((np.isfinite(sd))&(sd>-999)&(sd<9999)&(obs>-999)&(obs<9999))
                        err=sd[good]-obs[good]
                        mult=1
                        
                        RMS[i*len(scales)+j,t]=np.sqrt(np.mean((err)**2))
                        bias[i*len(scales)+j,t]=np.mean(err)
                        fullobs[i*len(scales)+j,t]=np.mean(obs[good])*mult
                        fullobs2[i*len(scales)+j,t]=np.mean(obs2[good])*mult
                    
                        print(" ".join([methods[i]+SDbase,scales[j],times[t],obj,
                                    str(RMS[i*len(scales)+j,t]),str(bias[i*len(scales)+j,t])]))
                    except Exception as e:
                        print(e)
    
        bias+=fullobs
        if obj=="dryspell":
            minmax=[0,6]
            # plot_vscale(bias,"dryspell",[-6,4],legend_ncol=2)
            # plot_vseason(bias,"dryspell",[-6,4],legend_ncol=2)
            if experiment=="conus":
                plot_vscale(bias,"dryspell",[0,8],legend_ncol=1,legend_loc=1,obs=fullobs,obs2=fullobs2,ylabel="Days")
                plot_vseason(bias,"dryspell",[0,15],legend_ncol=2,obs=fullobs,obs2=fullobs2,ylabel="Days")
            else:
                if pr_threshold==1:
                    plot_vscale(bias,"dryspell",[5,12],legend_ncol=1,legend_loc=1,obs=fullobs,obs2=fullobs2,ylabel="Days")
                    plot_vseason(bias,"dryspell",[2,25],legend_ncol=2,obs=fullobs,obs2=fullobs2,ylabel="Days")
                else:
                    plot_vscale(bias,"dryspell",[0,8],legend_ncol=1,legend_loc=1,obs=fullobs,obs2=fullobs2,ylabel="Days")
                    plot_vseason(bias,"dryspell",[0,15],legend_ncol=2,obs=fullobs,obs2=fullobs2,ylabel="Days")
            
        if obj=="wetspell":
            minmax=[1,15]
            # plot_vscale(bias,"wetspell",[-15,50])
            # plot_vseason(bias,"wetspell",[-15,50])
            if pr_threshold==1:
                plot_vscale(bias,"wetspell",[0,5],obs=fullobs,obs2=fullobs2,ylabel="Days")
                plot_vseason(bias,"wetspell",[0,10],obs=fullobs,obs2=fullobs2,ylabel="Days")
            else:
                plot_vscale(bias,"wetspell",[0,50],obs=fullobs,obs2=fullobs2,ylabel="Days")
                plot_vseason(bias,"wetspell",[0,40],obs=fullobs,obs2=fullobs2,ylabel="Days")
                
        if obj=="wetfrac":
            minmax=[0.05,0.6]
            # plot_vscale(bias,"wetfrac",[-0.15,0.7],legend_loc=1)
            # plot_vseason(bias,"wetfrac",[-0.15,1.0],legend_ncol=2)
            if pr_threshold==1:
                plot_vscale(bias,"wetfrac",[0.1,0.6],obs=fullobs,obs2=fullobs2,legend_loc=4,ylabel="")
                plot_vseason(bias,"wetfrac",[0.0,0.8],obs=fullobs,obs2=fullobs2,legend_ncol=2,ylabel="")
            else:
                plot_vscale(bias,"wetfrac",[0.1,1.0],obs=fullobs,obs2=fullobs2,legend_loc=4,ylabel="")
                plot_vseason(bias,"wetfrac",[0.2,1.3],obs=fullobs,obs2=fullobs2,legend_ncol=2,ylabel="")
                
        if obj=="MAP":
            minmax=[20,150]
            # plot_vscale(bias,"MAP",[-125,175],legend_ncol=3)
            # plot_vseason(bias,"MAP",[-20,50],multiplier=4.,legend_ncol=2)
            # plot_vseason(bias,"MAP",[0,100],obs=fullobs,obs2=fullobs2,multiplier=4.,legend_ncol=2,ylabel="mm")
            if experiment=="conus":
                plot_vscale(bias,"MAP",[500,1200],obs=fullobs,obs2=fullobs2,legend_ncol=3,ylabel="mm")
                plot_vseason(bias,"MAP",[0,150],obs=fullobs,obs2=fullobs2,legend_ncol=2,ylabel="mm")
            else:
                plot_vscale(bias,"MAP",[250,600],obs=fullobs,obs2=fullobs2,legend_ncol=3,ylabel="mm")
                plot_vseason(bias,"MAP",[0,100],obs=fullobs,obs2=fullobs2,legend_ncol=2,ylabel="mm")
        if obj=="interannual":
            minmax=[10,50]
            # plot_vscale(bias,"interannual",[-20,50],legend_ncol=2)
            # plot_vseason(bias,"interannual",[-20,30],multiplier=2.,legend_ncol=2)
            # plot_vseason(bias,"interannual",[0,50],obs=fullobs,obs2=fullobs2,multiplier=2.,legend_ncol=2,ylabel="mm")
            if experiment=="conus":
                plot_vscale(bias,"interannual",[50,250],obs=fullobs,obs2=fullobs2,legend_ncol=2,ylabel="mm")
                plot_vseason(bias,"interannual",[0,80],obs=fullobs,obs2=fullobs2,legend_ncol=2,ylabel="mm")
            else:
                plot_vscale(bias,"interannual",[50,150],obs=fullobs,obs2=fullobs2,legend_ncol=2,ylabel="mm")
                plot_vseason(bias,"interannual",[0,50],obs=fullobs,obs2=fullobs2,legend_ncol=2,ylabel="mm")
                
        bias-=fullobs
        make_plot(RMS,minmax,"RMS",obj,cmap=cm.Reds)
        minmax[0]=-minmax[1]
        make_plot(bias,minmax,"bias",obj,cmap=cm.seismic)
    
    yrs=["2yr","10yr","50yr","100yr"]
    # make extreme event plots
    nx=len(extremes)*4
    ny=len(methods)*len(scales)
    RMS=np.zeros((ny,nx))
    bias=np.zeros((ny,nx))
    fullobs=np.zeros((ny,nx))
    for o,obj in enumerate(extremes):
        for j in range(len(scales)):
            obs=io.read_nc("_".join([obs_base,scales[j],"annual",obj])+".nc").data
            if j==0:
                maskit=~np.isfinite(obs)
                obs=np.ma.array(obs,mask=maskit)
            if (o+j)==0:
                if plotmaps:
                    map_vis(obs[1,...],[obs_base,scales[j],"annual",obj],vmin=0,vmax=200)
            for i in range(len(methods)):
                try:
                    sd=io.read_nc("_".join([methods[i]+SDbase,scales[j],"annual",obj])+".nc").data
                    if j==0:
                        sd=np.ma.array(sd,mask=maskit|(~np.isfinite(sd)))
                    for curdim in range(len(sd.shape)):
                        if sd.shape[curdim]!=obs.shape[curdim]:
                            print("WARNING ERROR with shape for:"+ "_".join([methods[i]+SDbase,scales[j],"annual",obj]))
                            if curdim==0:
                                if obs.shape[curdim]<sd.shape[curdim]:
                                    sd=sd[:obs.shape[curdim],:]
                                else:
                                    obs=obs[:sd.shape[curdim],:]
                            if curdim==1:
                                if obs.shape[curdim]<sd.shape[curdim]:
                                    sd=sd[:,:obs.shape[curdim]]
                                else:
                                    obs=obs[:,:sd.shape[curdim]]
                            if curdim==2:
                                if obs.shape[curdim]<sd.shape[curdim]:
                                    sd=sd[:,:,:obs.shape[curdim]]
                                else:
                                    obs=obs[:,:,:sd.shape[curdim]]
                    if (o+j)==0:
                        if plotmaps:
                            map_vis(sd[1,...],[methods[i]+SDbase,scales[j],"annual",obj],vmin=0,vmax=200)
                            map_vis(sd[1,...]-obs[1,...],['error',methods[i]+SDbase,scales[j],"annual",obj],
                                    cmap=cm.seismic,vmin=-100,vmax=100)
                    for x in range(4):
                        err=sd[x,...]-obs[x,...]
                        RMS[i*len(scales)+j,o*4+x]=np.sqrt(np.mean(err**2))/np.mean(obs[x,...])
                        bias[i*len(scales)+j,o*4+x]=np.mean(err)/np.mean(obs[x,...])
                        fullobs[i*len(scales)+j,o*4+x]=np.mean(obs[x,...])
                        print(" ".join([methods[i]+SDbase,scales[j],obj,yrs[x],
                                str(RMS[i*len(scales)+j,o*4+x]),str(bias[i*len(scales)+j,o*4+x])]))
                except Exception as e:
                    print(e)

    minmax=[0.1,0.5]
    extreme_plot(RMS,minmax,"RMS",cmap=cm.Reds)
    minmax[0]=-minmax[1]
    extreme_plot(bias,minmax,"bias",cmap=cm.seismic)
    # plot_vscale(bias*100,"Percent Error Extremes",[-100,100],plotcol=2)
    if experiment=="conus":
        plot_vscale(bias*fullobs+fullobs,"Extremes",[0,200],legend_loc=1,plotcol=2,obs=fullobs,ylabel="mm/day 50yr return")
    else:
        plot_vscale(bias*fullobs+fullobs,"Extremes",[0,100],legend_loc=1,plotcol=2,obs=fullobs,ylabel="mm/day 50yr return")



