#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from stat_down import io
from mpl_toolkits.basemap import Basemap

# extremes=["extremes_nday1","extremes_nday2","extremes_nday3",
#           "extremes_nday4","extremes_nday5"]
objectives=["growing_season","frostdays","interannual_dtr","mean_dtr",
            "interannual_tmin","mean_tmin","interannual_tmax","mean_tmax",
            "interannual_tave","mean_tave"]
methods=["CAe0","CA","SARe0","SDe0","SDmon"]
# methods=["CAe0","CA","SDe0","SDmon"]
scales=["full_res","huc8","huc4","huc2"]
times=["annual","month01","month02","month03","month04","month05",
       "month06","month07","month08","month09","month10","month11",
       "month12","season1","season2","season3","season4"]

obs_base="obs-maurer.125-tasmax"
SDbase="-ncep-tasmax-BC12km"
methodnames=dict(CAe0="BCCAr",CA="BCCA",SARe0="SAR",SDe0="BCSDd",SDmon="BCSDm")

def map_vis(data,title=[],vmin=None,vmax=None,barlabel=None,outputdir=None):
    # print("visualizing : "+title)
    if len(data.shape)>2:
        plotdata=data[0,:,:]
    else:
        plotdata=data
    plt.clf()
    ny,nx=plotdata.shape
    geo=[35,43,-113,-101]
    m = Basemap(projection='merc',llcrnrlat=geo[0],urcrnrlat=geo[1],\
                llcrnrlon=geo[2],urcrnrlon=geo[3],lat_ts=(geo[0]+geo[1])/2)
    
    mapimg=m.imshow(plotdata,vmin=vmin,vmax=vmax)
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

def plot_vscale(data,plotname,yrange=None,plotcol=0,legend_loc=2,legend_ncol=1,obs=None,ylabel="Bias"):
    sdm=4;sdd=3;sar=2;cac=1;car=0
    varpos=[sdm,sdd,sar,cac,car]
    names=["BCSDm","BCSDd","SAR","BCCA","BCCAr"]
    colors=["red","darkred","green","blue","skyblue"]
    nscales=4
    
    plt.clf();
    ax=plt.gca()
    for v,n,c in zip(varpos,names,colors):
        ax.plot(data[v*nscales:(v+1)*nscales,plotcol],label=n,color=c,linewidth=2)
    if obs!=None:
        ax.plot(obs[0:nscales,plotcol],'-o',label="obs",color="k",linewidth=3)
    
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
    plt.savefig(plotname.replace(" ","_")+"_comp_scale.png")
    
def plot_vseason(data,plotname,yrange=None,plotrow=0,multiplier=1.0,legend_ncol=1,obs=None,ylabel="Bias"):
    sdm=4;sdd=3;sar=2;cac=1;car=0
    varpos=[sdm,sdd,sar,cac,car]
    names=["BCSDm","BCSDd","SAR","BCCA","BCCAr"]
    # colors=["red","magenta","green","blue","cyan"]
    colors=["red","darkred","green","blue","skyblue"]
    nscales=4
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
    plt.savefig(obj+"_"+name+".png")


    
if __name__ == '__main__':
    
    nx=len(times)
    ny=len(methods)*len(scales)
    
    # calculating maps for each statistics (objective function)
    for obj in objectives:
        RMS=np.zeros((ny,nx))
        bias=np.zeros((ny,nx))
        fullobs=np.zeros((ny,nx))
        for j in range(len(scales)):
            for t in range(len(times)):
                if ((obj=="frostdays") or (obj == "growing_season")) and times[t]!="annual":
                    pass
                else:
                    obs=io.read_nc("_".join([obs_base,scales[j],times[t],obj])+".nc").data
                    # if j==0:
                    #     map_vis(obs,[obs_base,scales[j],times[t],obj])
                    for i in range(len(methods)):
                        sd=io.read_nc("_".join([methods[i]+SDbase,scales[j],times[t],obj])+".nc").data
                        # if j==0:
                            # map_vis(sd,[methods[i]+SDbase,scales[j],times[t],obj])
                            # map_vis(sd-obs,["error",methods[i]+SDbase,scales[j],times[t],obj])
                        good=np.where((sd>-999)&(sd<9999)&(obs>-999)&(obs<9999))
                        err=sd[good]-obs[good]
                    
                        RMS[i*len(scales)+j,t]=np.sqrt(np.mean((err)**2))
                        bias[i*len(scales)+j,t]=np.mean(err)
                        fullobs[i*len(scales)+j,t]=np.mean(obs[good])
                        print(" ".join([methods[i]+SDbase,scales[j],times[t],obj,
                                    str(RMS[i*len(scales)+j,t]),str(bias[i*len(scales)+j,t])]))
    
        bias+=fullobs
        if obj[:len("mean")]=="mean":
            minmax=[0,2]
            if obj.split("_")[1]=="dtr":
                # plot_vscale(bias,obj,[15,20],obs=fullobs,legend_ncol=2,ylabel="K")
                plot_vseason(bias,obj,[13,24],obs=fullobs,legend_ncol=2,ylabel="K")
            if obj.split("_")[1]=="tmax":
                # plot_vscale(bias,obj,[16,20],obs=fullobs,legend_ncol=2,ylabel="K")
                # plot_vseason(bias,obj,[0,40],obs=fullobs,legend_ncol=2,ylabel="K")
                bias-=fullobs
                plot_vseason(bias,obj,[-2,2.5],legend_ncol=2,ylabel="bias [K]")
                bias+=fullobs
            if obj.split("_")[1]=="tave":
                # plot_vscale(bias,obj,[8,10],obs=fullobs,legend_ncol=2,ylabel="K")
                # plot_vseason(bias,obj,[-5,30],obs=fullobs,legend_ncol=2,ylabel="K")
                bias-=fullobs
                plot_vseason(bias,obj,[-2,2.5],legend_ncol=2,ylabel="bias [K]")
                bias+=fullobs
            if obj.split("_")[1]=="tmin":
                # plot_vscale(bias,obj,[-0.2,1.5],obs=fullobs,legend_ncol=2,ylabel="K")
                # plot_vseason(bias,obj,[-15,20],obs=fullobs,legend_ncol=2,ylabel="K")
                bias-=fullobs
                plot_vseason(bias,obj,[-2,2.5],legend_ncol=2,ylabel="bias [K]")
                bias+=fullobs
        if obj[:len("inter")]=="inter":
            minmax=[0,1]
            plot_vscale(bias,obj,[0.2,1.0],obs=fullobs,legend_ncol=2,ylabel="K")
            plot_vseason(bias,obj,[0,4],obs=fullobs,legend_ncol=2,ylabel="K")
        if obj=="growing_season":
            minmax=[10,30]
            plot_vscale(bias,"growing season",[70,140],obs=fullobs,legend_ncol=2,ylabel="days")
        if obj=="frostdays":
            minmax=[0,30]
            plot_vscale(bias,"frost days",[150,220],obs=fullobs,legend_ncol=2,ylabel="days")
        bias-=fullobs

        make_plot(RMS,minmax,"RMS",obj,cmap=cm.Reds)
        minmax[0]=-minmax[1]
        make_plot(bias,minmax,"bias",obj,cmap=cm.seismic)
        
    plot_histograms()
    


