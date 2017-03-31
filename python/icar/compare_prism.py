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
    
    icar_geo = mygis.read_geo(files[0])
    
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

def load_prism(dates):
    prism_data = []
    
    for d in dates:
        prism_date = d - datetime.timedelta(10)
        
        prism_file = master_prism.format(prism_date)
        prism_data.append( prism.load(prism_file) )
    
    return prism_data

def plot_comparison(case_name, icar, prism, date, geo, m=None):
    
    # bounds=[geo.lat.min(), geo.lat.max(), geo.lon.min(), geo.lon.max()]
    bounds=[37.7,38.25,-119.8,-119.1]
    
    datestring = str(date-datetime.timedelta(10))[:7]
    print("plotting "+datestring)
    
    plt.clf()
    plt.subplot(1,2,1)
    m = vis(icar, lat=geo.lat, lon=geo.lon, reproject=True, geo=bounds, 
            title="ICAR "+datestring, latstep=0.25, lonstep=0.25, m=m,clim=(0,800))
    
    plt.subplot(1,2,2)
    m = vis(prism.data, lat=prism.geo.lat, lon=prism.geo.lon, reproject=True, geo=bounds, 
            title="PRISM "+datestring, latstep=0.25, lonstep=0.25, m=m,clim=(0,800))#clim=(icar.min(), icar.max()))
    
    plt.savefig(case_name+"_"+datestring+".png")
    
    return m

def main(test_case):
    icar_data  = load_data(test_case)
    prism_data = load_prism(icar_data.dates)
    
    plt.figure(figsize=(13,4))
    m=None
    for i in range(len(prism_data)):
        m = plot_comparison(test_case, icar_data.data[i], prism_data[i], icar_data.dates[i], icar_data.geo, m=m)

if __name__ == '__main__':
    if len(sys.argv)<2:
        print("You must specify an input case (directory name)")
    else:
        main(sys.argv[1])
