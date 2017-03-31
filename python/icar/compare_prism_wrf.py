#!/bin/env python
import sys
from glob import glob
import datetime

import numpy as np
import matplotlib.pyplot as plt

import prism
from bunch import Bunch
import mygis
from stat_down.map_vis import vis

master_prism = "/glade/scratch/gutmann/prism/PRISM_ppt_stable_4kmM3_{0.year}{0.month:02}_bil.bil"
master_wrf = "/glade/p/ral/RHAP/asd000/CONUS/CTRL_new/monthly/2D/wrf2d_total/wrf2d_d01_monthly_total_{0.year}{0.month:02}.nc"
wrf_conus_geo_file = "/glade/p/ral/RHAP/asd000/CONUS/wrfout_conus_constants.nc"
geo_file = "subset_geo_em.nc"

def read_date_from_file(filename):
    dates = filename.split("_")[-4:]
    dates[-1] = dates[-1][:2] # drop the potential ".nc" at the end
    dates = np.array(dates,dtype=int)
    
    return datetime.datetime(*dates)

def load_data(test_case):
    files = glob(test_case+"/output/short_icar_out_20??_??_01*")
    files.sort()
    files2=glob(test_case+"/output/icar_out_20??_??_01*")
    files2.sort()
    
    files.extend(files2)
    
    icar_geo = mygis.read_geo(test_case+"/"+geo_file)
    
    data = []
    last_data = mygis.read_nc(files[0],"rain").data[0]
    for f in files[1:]:
        this_data = mygis.read_nc(f,"rain").data[0]
        data.append( this_data - last_data )
        last_data = this_data
    
    dates = []
    for f in files[1:]:
        dates.append( read_date_from_file(f) )

    return Bunch(data=data, dates=dates, geo=icar_geo)

def load_wrf(dates):
    wrf_data=[]
    
    wrf_geo = mygis.read_geo(wrf_conus_geo_file)
    
    for d in dates:
        try:
            wrf_date = d - datetime.timedelta(10)
            
            wrf_file = master_wrf.format(wrf_date)
            curdata = mygis.read_nc(wrf_file,"PREC_ACC_NC").data[0]
            
            wrf_data.append(Bunch(data=curdata, geo=wrf_geo) )
        except:
            print("Could not load wrf for date: "+str(d))
            wrf_data.append(Bunch(data=None, geo=wrf_geo) )
    
    return wrf_data

def load_prism(dates):
    prism_data = []
    
    for d in dates:
        try:
            prism_date = d - datetime.timedelta(10)
            
            prism_file = master_prism.format(prism_date)
            prism_data.append( prism.load(prism_file) )
        except:
            print("Could not load PRISM for date: "+str(d))
    
    return prism_data

def plot_comparison(case_name, icar, prism, wrf, date, geo, m=None):
    
    bounds=[geo.lat.min(), geo.lat.max(), geo.lon.min(), geo.lon.max()]
    # bounds=[37.7,38.25,-119.8,-119.1]
    # bounds[1]=52
    # bounds[2]=-125
    # bounds[3]=-108
    
    datestring = str(date-datetime.timedelta(10))[:7]
    print("plotting "+datestring)
    
    latstep = 0.25
    lonstep = 0.25
    
    plt.clf()
    
    plt.subplot(1,3,1)
    m = vis(prism.data, lat=prism.geo.lat, lon=prism.geo.lon, reproject=True, geo=bounds, colorbar=False,
            title="PRISM "+datestring, latstep=latstep, lonstep=lonstep, m=m,clim=(0,500))#clim=(icar.min(), icar.max()))

    plt.subplot(1,3,2)
    m = vis(icar, lat=geo.lat, lon=geo.lon, reproject=True, geo=bounds, colorbar=False,
            title="ICAR "+datestring, latstep=latstep, lonstep=lonstep, m=m,clim=(0,500))
            
    plt.subplot(1,3,3)
    m = vis(wrf.data, lat=wrf.geo.lat, lon=wrf.geo.lon, reproject=True, geo=bounds, colorbar=False,
            title="WRF "+datestring, latstep=latstep, lonstep=lonstep, m=m,clim=(0,500))#clim=(icar.min(), icar.max()))
    
    plt.savefig("prism_wrf_icar_"+case_name+"_"+datestring+".png")
    
    return m

def main(test_case):
    icar_data  = load_data(test_case)
    prism_data = load_prism(icar_data.dates)
    wrf_data = load_wrf(icar_data.dates)
    
    plt.figure(figsize=(18,4))
    m=None
    for i in range(len(wrf_data)):
        m = plot_comparison(test_case, icar_data.data[i], prism_data[i], wrf_data[i], icar_data.dates[i], icar_data.geo, m=m)

if __name__ == '__main__':
    if len(sys.argv)<2:
        print("You must specify an input case (directory name)")
    else:
        main(sys.argv[1])
