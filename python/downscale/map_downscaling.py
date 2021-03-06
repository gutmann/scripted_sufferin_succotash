#!/usr/bin/env python
import numpy as np
import glob
from bunch import Bunch
import mygis

cur_file="prob_stn*_prcp_norm.nc"
global_output_file="grid_"+"prob_prcp_norm.nc"
# cur_file="stn*_prcp_norm.nc"
# global_output_file="grid_"+"prcp_norm.nc"
global_basefile="/d2/gutmann/usbr/stat_data/DAILY/obs/maurer.125/pr/nldas_met_update.obs.daily.pr.2000.nc"
# data_dir="processing_*/"
data_dir=""
global_file_search=data_dir+cur_file
global_bounds=Bunch(lat=Bunch(min=34,max=44),lon=Bunch(min=-114,max=-100))


def find_files(startdate):
    """find the files that have "run_start_date" equal to startdate"""
    all_files=glob.glob(global_file_search)
    output_files=[]
    for f in all_files:
        date_data=mygis.read_nc(f,"run_start_date").data
        thisdate="".join(date_data[0].data)
        if thisdate==startdate:
            output_files.append(f)
    
    print(output_files)        
    return output_files

def load_geo(filename):
    """load a 2d lat and lon grid from filename"""
    latnames=["lat","latitude","XLAT"]
    lonnames=["lon","longitude","XLONG"]
    lat=None
    lon=None
    
    for l in latnames:
        if lat==None:
            try:
                lat=mygis.read_nc(filename,l).data
            except:
                lat=None
    for l in lonnames:
        if lon==None:
            try:
                lon=mygis.read_nc(filename,l).data
            except:
                lon=None

    ymin=np.where(lat>global_bounds.lat.min)[0][0]
    ymax=np.where(lat>global_bounds.lat.max)[0][0]
    lat=lat[ymin:ymax]
    
    xmin=np.where(lon>global_bounds.lon.min)[0][0]
    xmax=np.where(lon>global_bounds.lon.max)[0][0]
    lon=lon[xmin:xmax]

    if len(lon.shape)==1:
        lon,lat=np.meshgrid(lon,lat)
        
    return Bunch(lon=lon,lat=lat)
        



def get_xy(filename,geo):
    """return x,y coordinates in geo for all points in the file"""
    lat=mygis.read_nc(filename,"station_latitude").data
    lon=mygis.read_nc(filename,"station_longitude").data
    xy=Bunch(x=np.zeros(len(lon)),y=np.zeros(len(lat)))
    for i in range(lat.size):
        dists=(lat[i]-geo.lat)**2+(lon[i]-geo.lon)**2
        y,x=np.unravel_index(np.argmin(dists),dists.shape)
        xy.x[i]=x
        xy.y[i]=y
        
    return xy
    

def main(start_date="19900101"):
    """convert a file (or files) from downscaling code to maps"""
    files=find_files(start_date)
    files.sort()
    geo=load_geo(global_basefile)
    times=mygis.read_nc(files[0],"time").data
    ntimes=times.size
    
    tmp=mygis.read_nc(files[0],"coefficient",returnNCvar=True)
    output_data=np.zeros((tmp.data.shape[1],ntimes,geo.lon.shape[0],geo.lon.shape[1]))
    tmp.ncfile.close()
    print(output_data.shape)
    for f in files:
        print(f)
        data=mygis.read_nc(f,"coefficient").data
        locations=get_xy(f,geo)
        for i in range(len(locations.x)):
            output_data[:,:,locations.y[i],locations.x[i]]=data[0,:,:,i]
        
    print("Writing output file")
    mygis.write(global_output_file,output_data,varname="coefficient",dtype="d",dims=("variable","time","latitude","longitude"),
                    extravars=[ Bunch(data=times,name="time",dims=("time",),dtype="d",
                                      attributes=Bunch(units="seconds since 1970-01-01 00:00:00.0 0:00")),
                                Bunch(data=geo.lat,name="latitude",dims=("latitude","longitude"),dtype="f",
                                      attributes=Bunch(units="degrees")),
                                Bunch(data=geo.lon,name="longitude",dims=("latitude","longitude"),dtype="f",
                                      attributes=Bunch(units="degrees"))])

    
        


if __name__ == '__main__':
    main()