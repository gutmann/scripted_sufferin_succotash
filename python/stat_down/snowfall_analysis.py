#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import mygis

start_current_year = 1995
end_current_year   = 2005
start_future_year  = 2045
end_future_year    = 2055

max_temp_threshold =  2 # upper snowfall bound, temperatures above this value have 0% snow
min_temp_threshold = -2 # lower rainfall bound, temperatures below this value have 0% rain

ofile       = "month_{1:02}_{0}_{2}_snowfall.nc"
ifile       = "forcing/{0}/{0}_Forcing.{1}-{2:02}.nc"
plotfile    = "snowplot_{0}_month_{1:02}.png"
delta_title = "Fut-Cur {0} month:{1} {2}"

def compute_mean_snowfall(method,startyear,endyear,month):
    """docstring for compute_snowfall"""
    snow_total = None
    prec_total = None
    
    n = (endyear-startyear)+1
    for y in range(startyear,endyear+1):
        # print(ifile.format(method,y,month))
        precip = mygis.read_nc(ifile.format(method,y,month),"ppt").data  * 3600.0 # convert to mm/hr (from mm/s)
        temp   = mygis.read_nc(ifile.format(method,y,month),"temp").data - 273.15 # convert to deg C
        
        snow_ratio = (temp - min_temp_threshold)
        if (max_temp_threshold-min_temp_threshold)>0:
            snow_ratio /= (max_temp_threshold-min_temp_threshold)
        else:
        snow_ratio=np.clip(1-snow_ratio,0,1)
        
        snow=precip*snow_ratio
        if snow_total==None:
            snow_total  = snow.sum(axis=0)
            prec_total  = precip.sum(axis=0)
        else:
            snow_total += snow.sum(axis=0)
            prec_total += precip.sum(axis=0)
            
    
    snow_total /= n
    prec_total /= n
    return snow_total,prec_total


def plot_data(current,future,precip_current,precip_future, method,month):
    """docstring for plot_data"""
    plt.figure(figsize=(15,8))
    
    plt.subplot(231)
    plt.imshow(future-current,cmap=plt.cm.seismic_r)
    plt.clim(-100,100)
    plt.title(delta_title.format(method,month,"Snow"))
    plt.colorbar()

    plt.subplot(232)
    plt.imshow(current)
    plt.title("Current Snowfall")
    plt.clim(0,300)
    plt.colorbar()

    plt.subplot(233)
    plt.imshow(future)
    plt.title("Future Snowfall")
    plt.clim(0,300)
    plt.colorbar()
    
    plt.subplot(234)
    plt.imshow(precip_future-precip_current,cmap=plt.cm.seismic_r)
    plt.clim(-100,100)
    plt.title(delta_title.format(method,month,"Precip"))
    plt.colorbar()

    plt.subplot(235)
    plt.imshow(precip_current)
    plt.title("Current Precip")
    plt.clim(0,300)
    plt.colorbar()

    plt.subplot(236)
    plt.imshow(precip_future)
    plt.title("Future Precip")
    plt.clim(0,300)
    plt.colorbar()

    plt.savefig(plotfile.format(method,month))
    plt.close()
    

def main():
    """Loop through Downscaling Methods and Months computing monthly snowfall, writing to disk and making summary plots"""
    # for method in ["BCCA12K","BCSAR12K","BCSD12K","BCSDdisag12K"]:
    for method in ["BCSAR12K","BCSD12K","BCSDdisag12K"]:
        for month in range(1,13):
            print(method, month)
            csnow,cprec=compute_mean_snowfall(method, start_current_year, end_current_year, month)
            fsnow,fprec=compute_mean_snowfall(method, start_future_year, end_future_year, month)
        
            # mygis.write(ofile.format(method,month,"current"),current,varname="snowfall")
            # mygis.write(ofile.format(method,month,"future"),future,varname="snowfall")
        
            plot_data(csnow,fsnow,cprec,fprec,method,month)
            # except:
            #     print("Error working with : "+method+" "+str(month))
            
        

if __name__ == '__main__':
    main()