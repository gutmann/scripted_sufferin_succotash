#!/usr/local/bin/python
import matplotlib.pyplot as plt
import numpy as np
import glob
import mygis
# import map
model_names={"cnrm":"CNRM-CM5","ccsm":"CCSM4","cesm":"CCSM4","gfdl":"GFDL-CM3","miro":"MIROC5","nore":"NorESM"}
model_order={"cnrm":1,"ccsm":2,"gfdl":3,"miro":4,"nore":5}

base_file="daily_data/icar_run_era/icar_1990_daily_rain.nc"
era_annual_p="daily_data/icar_run_era/annual_mean_precip.nc"
cmip_hist = "CMIP5/summary_*_historical_precip.nc.nc"

rcp45_file="daily_data/icar_run_*_rcp85/annual_mean_precip.nc"
cmip_rcp45 = "CMIP5/summary_*_rcp85_precip.nc.nc"

def main():
    geo=mygis.read_geo(base_file)
    erap=mygis.read_nc(era_annual_p).data
    rcp45_files=glob.glob(rcp45_file)
    print(rcp45_file, rcp45_files)
    rcp45_files.sort()
    iens_h=mygis.read_files(rcp45_files)
    imean=np.zeros(iens_h[0].shape)
    
    cmip_hist_files=glob.glob(cmip_hist)
    cmip_hist_files.sort()
    pens_h=mygis.read_files(cmip_hist_files,"precip")
    pmean=np.zeros(pens_h[0].shape)
    for p in pens_h:
        pmean+=p
    pmean/=len(pens_h)
    
    pens_r_files=glob.glob(cmip_rcp45)
    pens_r_files.sort()
    pens_r=mygis.read_files(pens_r_files,"precip")
    plt.figure(figsize=(8,4))
    m=None
    for i in range(len(pens_r_files)):
        model=pens_r_files[i].split("/")[-1].split("_")[1][:4].lower()
        print(model)
        
        plt.clf()
        plt.subplot(1,2,1)
        print(pens_r[i].shape, pmean.shape)
        # m=map.vis(pens_r[i][-1] - pmean[-1],lat=geo.lat, lon=geo.lon, 
        #          cmap=plt.cm.seismic_r, m=m)
        plt.imshow(pens_r[i][-1] - pmean[-1],cmap=plt.cm.seismic_r)
        plt.clim(-500,500)
        plt.title(model_names[model])
        model1=model
        
        model=rcp45_files[i].split("/")[-2].split("_")[2][:4].lower()
        plt.subplot(1,2,2)
        print(iens_h[i].shape)
        plt.imshow(iens_h[i]-erap,cmap=plt.cm.seismic_r)
        # m=map.vis(iens_h[i]-erap,lat=geo.lat, lon=geo.lon, 
        #         cmap=plt.cm.seismic_r, m=m)
        plt.clim(-500,500)
        plt.title(model_names[model])
        # plt.tight_layout()
        plt.draw()
        plt.savefig("{}_{}_deltas_rcp85.png".format(model1.upper(), model.upper()))
        

if __name__ == '__main__':
    main()
