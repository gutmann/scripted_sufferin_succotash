#!/usr/bin/env python
import numpy as np
from stat_down.map_vis import vis as map_vis
import mygis
import matplotlib.pyplot as plt

# CCSM file start date=1960,enddate = 2080
# 20th century dates
c20start = 35*365 # jan1 1995
c20end   = 45*365 # jan1 2005
# 21st century dates
c21start = 85*365 # Jan1 2045
c21end   = 95*365 # Jan1 2055
title_20c=" Precip 1995-05 (mm)"
title_21c=" Precip 2045-55 (mm)"

# Longer averaging time
# 20th century dates
c20start = 15*365 # jan1 1975
c20end   = 45*365 # jan1 2005
# 21st century dates
c21start = 75*365 # Jan1 2035
c21end   = 105*365 # Jan1 2065
title_20c=" Precip 1975-05 (mm)"
title_21c=" Precip 2035-65 (mm)"

precipmax=30
# precipmax=150
rbound=1.5

xlim=(235,295)
ylim=(25,53)

monthdays=[31,28,31,30,31,30,31,31,30,31,30,31]
monthnames=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
monthstarts=np.zeros(12)
monthends=np.zeros(12)
for i in range(12):
    if i==0:
        monthstarts[i]=0
    else:
        monthstarts[i]=monthends[i-1]
    monthends[i]=monthstarts[i]+monthdays[i]
    
ccsm_file="prate.run5.complete.1960-2080.nc"
def main():
    """visualize each month of ccsm precip data"""
    d=mygis.read_nc(ccsm_file,"pr",returnNCvar=True)
    c20=d.data[c20start:c20end,...]
    c21=d.data[c21start:c21end,...]
    d.ncfile.close()
    
    geo=mygis.read_geo(ccsm_file)
    geo.lon[geo.lon<0]+=360 #convert range [-180 - 0] into [180-360]
    
    doy=np.arange(c20.shape[0]) % 365
    plt.figure(figsize=(4*1.5,7*1.5))
    for month in range(12):
        print(monthnames[month])
        cur=np.where((doy>monthstarts[month]) & (doy<monthends[month]))
        cur20=c20[cur[0],...].mean(axis=0)*monthdays[month]*86400
        cur21=c21[cur[0],...].mean(axis=0)*monthdays[month]*86400
        
        plt.clf();
        plt.subplot(311)
        map_vis(cur20,vmin=0,vmax=precipmax,title=monthnames[month]+title_20c,
                geo=None,lat=geo.lat,lon=geo.lon,xlim=xlim,ylim=ylim,cmap=plt.cm.jet_r)
        plt.subplot(312)
        map_vis(cur21,vmin=0,vmax=precipmax,title=monthnames[month]+title_21c,
                geo=None,lat=geo.lat,lon=geo.lon,xlim=xlim,ylim=ylim,cmap=plt.cm.jet_r)

        plt.subplot(313)
        ratio=cur21/cur20
        lower=ratio<1
        ratio[lower]=0-(1/ratio[lower])
        ratio[(ratio>0) & (ratio< rbound)]=0
        ratio[(ratio<0) & (ratio>-rbound)]=0
        # map_vis(ratio,vmin=0,vmax=2,title=monthnames[month]+
        #         " Ratio [fut/cur]",
        map_vis(ratio,vmin=-5,vmax=5,title=monthnames[month]+
                " Ratio fut/cur (fut>cur) & 0-cur/fut (cur>fut)",
                geo=None,lat=geo.lat,lon=geo.lon,xlim=xlim,ylim=ylim,
                cmap=plt.cm.seismic_r)
        
        plt.savefig(str(month+1)+"_summary_"+monthnames[month]+".png")
    
    

if __name__ == '__main__':
    main()