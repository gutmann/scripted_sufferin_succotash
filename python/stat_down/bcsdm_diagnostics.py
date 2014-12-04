#!/usr/bin/env python
from __future__ import print_function
import os,sys
import datetime as dt

import matplotlib.pyplot as plt
import numpy as np

import mygis


cur_start_time = dt.datetime(1995,10,1)
cur_end_time   = dt.datetime(2006,10,1)
fut_start_time = dt.datetime(2045,10,1)
fut_end_time   = dt.datetime(2056,10,1)

month_lengths=np.array([31,28.25,31,30,31,30,31,31,30,31,30,31])
output_dir="/glade/u/home/gutmann/scratch/usbr/stat_data/cc_test/bcsdm_diagnosis/forcing/"
output_dir="./"
# data stored in : /glade/p/ral/RHAP/naoki/hydro_ccsm/vic/models/vic/output/netcdf/BCSDdisag12K/UCO
xlabels=["J","F","M","A","M","J","J","A","S","O","N","D"]

def plot_runoff_delta(varname="runoff",extra_var=None,aggregation=np.sum,monthly_total=True,multiplier=86400.0,method="",calendar="noleap"):
    """docstring for plot_runoff_delta"""
    print("Loading data :"+varname)
    # read in runoff from all of the netcdf files in the current directory and concatenate them along axis 0 (time)
    data  = mygis.read_files(method+"_*[1-2]*.nc",varname,axis=0)*multiplier # mm/s * 86400 s/day = mm/day)
    if (extra_var!=None):
        data += mygis.read_files(method+"_[1-2]*.nc",extra_var,axis=0)*multiplier # mm/s * 86400 s/day = mm/day)
    data=data[:,60:150,80:200]
    # create a list of datetime objects to match the netcdf data
    if calendar=="noleap":
        # to fake a 365 day calendar
        times=[dt.datetime(1980,10,1,0,0)+dt.timedelta(np.floor(i/365.0*365.25)) for i in range(data.shape[0])]
    else:
        times=[dt.datetime(1960,1,1,0,0)+dt.timedelta(i) for i in range(data.shape[0])]
    
    # summary data (12 months of no data to begin with)
    n_cur = [0]*12
    curdata=[None]*12
    n_fut = [0]*12
    futdata=[None]*12
    print("Computing Monthly Means")
    for month in range(1,13):
        print("  Month:{0:02}".format(month),end="\r")
        sys.stdout.flush()
        for i,t in enumerate(times):
            # if this time step is the month we are analyzing:
            if t.month==month:
                # if this time step is in the current time period:
                if (t>=cur_start_time) and (t<cur_end_time):
                    if (curdata[month-1]==None):
                        curdata[month-1] = data[i,:,:]
                        n_cur[month-1]   = 1
                    else:
                        curdata[month-1] += data[i,:,:]
                        n_cur[month-1]   += 1

                # if this time step is in the future time period:
                if (t>=fut_start_time) and (t<fut_end_time):
                    if (futdata[month-1]==None):
                        futdata[month-1] = data[i,:,:]
                        n_fut[month-1]   = 1
                    else:
                        futdata[month-1] += data[i,:,:]
                        n_fut[month-1]   += 1
    
    print("Plotting    ")
    monthly_deltas=np.zeros(12)
    for i in range(12):
        print("  Month:{0:02} n={1},{2}    ".format(i+1,n_cur[i],n_fut[i]),end="\r")
        sys.stdout.flush()
        # divide by n to get daily mean and multiply by days/month to get monthly totals
        curdata[i] /= n_cur[i]
        futdata[i] /= n_fut[i]
        if monthly_total:
            curdata[i]*=month_lengths[i]
            futdata[i]*=month_lengths[i]
        
        # make a map of the current monthly changes
        plt.clf()
        if varname[0].lower()=="t":
            cmap=plt.cm.seismic
        else:
            cmap=plt.cm.seismic_r    
        plt.imshow(futdata[i] - curdata[i],cmap=cmap)
        max_range=np.max(np.abs(futdata[i] - curdata[i]))
        plt.clim(-max_range,max_range)
        # plt.clim(-150,150)
        plt.colorbar()
        plt.title(varname+"_"+method+" delta month{0:02}".format(i+1))
        plt.savefig(output_dir+varname+"_"+method+"_delta_month{0:02}.png".format(i+1))
        
        # compute the sum of all runoff changes in the basin (the mask isn't necessary, was intended to mask out points outside the basin)
        # mask=np.where(futdata[i]!=0)
        monthly_deltas[i]=aggregation((futdata[i] - curdata[i]))  #[mask])
        
    # make a plot of the basin total changes
    plt.clf()
    plt.plot(monthly_deltas)
    # print(monthly_deltas)
    plt.plot([0,12],[0,0],color="black")
    plt.xticks(range(12),xlabels)
    plt.savefig(output_dir+varname+"_"+method+"_monthly_change.png")
    print("Finished              ")
    print("-"*20)
        
    
# def plot_time_series(varname,method,extra_var=None,monthly_total=True):
#     """docstring for plot_time_series"""
#     print("Loading data :"+varname)
#     # read in runoff from all of the netcdf files in the current directory and concatenate them along axis 0 (time)
#     data  = mygis.read_files(method+"_[1-2]*.nc",varname,axis=0)*multiplier # mm/s * 86400 s/day = mm/day)
#     if (extra_var!=None):
#         data += mygis.read_files(method+"_[1-2]*.nc",extra_var,axis=0)*multiplier # mm/s * 86400 s/day = mm/day)
#     data=data[:,60:150,80:200]
#     # create a list of datetime objects to match the netcdf data
#     if calendar=="noleap":
#         # to fake a 365 day calendar
#         times=[dt.datetime(1980,10,1,0,0)+dt.timedelta(np.floor(i/365.0*365.25)) for i in range(data.shape[0])]
#     else:
#         times=[dt.datetime(1979,1,1,0,0)+dt.timedelta(i) for i in range(data.shape[0])]
    
    


def main():
    """docstring for main"""
    # method=os.getcwd().split("/")[-2]
    # print(method)
    # plot_runoff_delta(varname="runoff",extra_var="baseflow",aggregation=np.mean,method=method)
    # plot_runoff_delta("swe",aggregation=np.mean,monthly_total=False,multiplier=1.0,method=method)
    # plot_runoff_delta("prec",aggregation=np.mean,method=method)
    
    method=os.getcwd().split("/")[-1]
    method="BCSDdisag_12km"
    # plot_runoff_delta("Prec",aggregation=np.mean,method=method,multiplier=1.0,calendar="gregorian")
    plot_runoff_delta("tasmax",aggregation=np.mean,method=method,multiplier=1.0,calendar="gregorian",monthly_total=False)
    # plot_runoff_delta("pr",aggregation=np.mean,method=method,multiplier=1.0,calendar="gregorian",monthly_total=False)
    # plot_runoff_delta("Tmax",aggregation=np.mean,method=method,multiplier=1.0,calendar="gregorian",monthly_total=False)
    # plot_runoff_delta("Tmin",aggregation=np.mean,method=method,multiplier=1.0,calendar="gregorian",monthly_total=False)
    
    # plot_time_series("Tmin",method=method,monthly_total=False)
    # plot_time_series("Tmax",method=method,monthly_total=False)
    # plot_time_series("Prec",method=method)

if __name__ == '__main__':
    main()

