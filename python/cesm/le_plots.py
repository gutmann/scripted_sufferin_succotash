import numpy as np
import matplotlib.pyplot as plt


def make_plots(plot_fun, data):
    """docstring for make_plots"""
    
    plt.figure(figsize=(30,35))
    for i in range(len(data)):
        plt.subplot(7,6,i+1)
        plot_fun(data[i])
        plt.title("Ens:{:03}".format(i+1),fontsize=18)
    plt.tight_layout()

def timeseries(data):
    """docstring for timeseries"""
    times=data.years
    ens_i=data.raw
    ens_m=data.ens_mean
    
    if "fill" in data.keys():
        for d in data.fill:
            plt.fill_between(d[0],0,d[1],color="black",alpha=0.25)
    
    print(times.shape,ens_m.shape)
    n=min(times.shape[0],ens_m.shape[0])
    plt.plot(times[:n],ens_m[:n],color="black", linewidth=2)
    n=min(times.shape[0],ens_i.shape[0])
    plt.plot(times[:n],ens_i[:n],color="black", linewidth=1,alpha=0.25)
    plt.xlim(data.daterange)
    plt.ylim(data.yrange)
    
    if "smooth" in data.keys():
        plt.plot(times,data.smooth,color="black",linewidth=1)

    if "fillpts" in data.keys():
        plt.plot(data.fillpts[0],data.fillpts[1],'D', markersize=8, color="black")

    if "cool_smooth" in data.keys():
        print("plotting cool_smooth")
        plt.plot(times,data.cool_smooth, color="blue")
    if "cool_mean" in data.keys():
        print("plotting cool_mean")
        plt.plot(times,data.cool_mean, color="blue",linewidth=2)
    if "coolpts" in data.keys():
        plt.plot(data.coolpts[0],data.coolpts[1],'D', markersize=8, color="blue")

    if "warm_smooth" in data.keys():
        print("plotting warm_smooth")
        plt.plot(times,data.warm_smooth, color="red")
    if "warm_mean" in data.keys():
        print("plotting warm_mean")
        plt.plot(times,data.warm_mean, color="red",linewidth=2)
    if "warmpts" in data.keys():
        plt.plot(data.warmpts[0],data.warmpts[1],'D', markersize=8, color="red")
    
    
    
def histogram(data):
    """docstring for histogram"""
    xvals=data[0]
    yvals=data[1]
    x0=data[2]
    # xlim=data[3]
    # ylim=data[4]
    for i in range(len(yvals)):
        # x=[xvals[i]-dx/2,xvals[i]-dx/2,xvals[i]+dx/2,xvals[i]+dx/2]
        # y=[0,yvals[i],yvals[i],0]
        x=[xvals[i],  xvals[i+1]]
        y=[yvals[i],  yvals[i]]
        plt.fill_between(x,0,y,color="grey")
    
    plt.plot([x0,x0],[0,yvals.max()+1],color="red",linewidth=3)

