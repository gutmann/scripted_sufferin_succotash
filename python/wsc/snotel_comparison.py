#!/usr/bin/env python
import glob,re
import numpy as np

from bunch import Bunch
import swim_io as ncio


snotel_list_filename="snotel/sitelist.txt"
wrf_geo_file="wrf/4km_wrf_output.nc"
wrf_data_file="wrf/wrf_prcp_TIMEERROR.nc"

def load_snotel_list(snotel_list_file):
    """Load station name/location information from a text file
    
    example file:
ID     Site Name         Latitude  Longitude   Elevation
05J04S PHANTOM VALLEY    40.38333 -105.83333     2752.4
05J05S WILD BASIN    40.20000 -105.60000     2913.9
05J06S DEADMAN HILL    40.80000 -105.76667     3115.1
    """
    with open(snotel_list_file) as f:
        for i in range(1):l=f.next()
        
        outputlist=[]
        for l in f:
            try:
                line=l.split()
                name=line[0].strip()
                lat=np.double(line[-3])
                lon=np.double(line[-2])
                outputlist.append(Bunch(name=name,lat=lat,lon=lon))
            except Exception as e: 
                print(e)
                print(l)
        
    return outputlist

def load_snotel_data(snotel_station):
    """load data from a snotel file 
    file format (originally Kyoko?)
    
name = BUTTE
index = 06L11S
lat. = 38.87250
lon. = -106.95278
elev. = 3096.8
wrf1 =   1 -106.9525   38.8705  2888.30  244  303
wrf2 =   2 -106.9525   38.8890  3072.10  245  303
wrf3 =   3 -106.9762   38.8705  2720.00  244  302
wrf4 =   4 -106.9287   38.8704  3018.00  244  304
JulianDay  Acc.Precip(mm)
     2444878    0.00
     2444879    0.00
    """
    try:
        filename=glob.glob("snotel/"+snotel_station.name+"_data.txt")[0]
        with open(filename) as f:
            l=""
            while (not re.match(".*JulianDay .*",l)):l=f.next()
        
            outputlist=[]
            lastprecip=0.0
            for l in f:
                line=l.split()
                jday=int(line[0])
                precip=float(line[1])-lastprecip
                if precip<-100:
                    precip=float(line[1])
                elif precip<0:
                    precip=0
                outputlist.append([jday,precip])
                lastprecip=float(line[1])
        return np.array(outputlist)
        
    except Exception as e:
        print("ERROR:"+str(e))
        print("  With station:"+snotel_station.name)
        return -1

def find_wrf_point(geo,point):
    """find array index in geo for point"""
    dists=(geo.lat-point.lat)**2+(geo.lon-point.lon)**2
    return np.unravel_index(np.argmin(dists),geo.lat.shape)

def load_wrf_data(station,wrfgeo,wrfdatafile,data=None):
    """load the time series of wrf precip data from wrfdatafile for the GHCN station given"""
    
    (i,j)=find_wrf_point(wrfgeo,station)
    if (data==None):
        wrfdata=ncio.read_nc(wrfdatafile,"prec",returnNCvar=True)
        outputdata=wrfdata.data[:,i,j] #its not clear if NIO doesn't support this, something is breaking here. 
        wrfdata.ncfile.close()
    else:
        outputdata=data[:,i,j]
    return outputdata


def load_wrf_geo(filename):
    """load WRF lat/lon data from netcdf file"""
    
    lat=ncio.read_nc(filename,"XLAT").data[0,...] #index 0 in the time dimension
    lon=ncio.read_nc(filename,"XLONG").data[0,...]
    return Bunch(lat=lat,lon=lon)

def calc_stats(data):
    data=data[np.where(data>=0)]
    mean_annual=data.mean()*365.25
    wetdays0=np.zeros(data.size)
    wetdays0[data>0]=1
    wetdays1=np.zeros(data.size)
    wetdays1[data>1]=1
    wetdays2=np.zeros(data.size)
    wetdays2[data>2.54]=1
    wetfrac0=wetdays0.sum()/float(data.size)
    wetfrac1=wetdays1.sum()/float(data.size)
    wetfrac2=wetdays2.sum()/float(data.size)
    
    data=data[data>0]
    if (data.size>100):
        data.sort()
        p99=data[int(round(data.size*0.99))]
        p95=data[int(round(data.size*0.95))]
        f99=data[int(round(data.size*0.99)):].sum()/data.sum()
        f95=data[int(round(data.size*0.95)):].sum()/data.sum()
    else:
        p99=-999
        p95=-999
        f99=-999
        f95=-999
    
    return Bunch(mean=mean_annual,wetfrac0=wetfrac0,wetfrac1=wetfrac1,wetfrac2=wetfrac2,p99=p99,p95=p95,f99=f99,f95=f95)

def write_stats(stats):
    """docstring for write_stats"""
    with open("snotel_stats.txt","w") as f:
        f.write("lat lon ")
        for k in stats[0].wrf.keys():
            f.write(k+" ")
        f.write("\n")
        for stat in stats:
            f.write(str(stat.station.lat)+" "+str(stat.station.lon))
            for k in stat.wrf.keys():
                f.write(" ")
                f.write(str(stat.snotel[k])+" "+str(stat.wrf[k]))
            f.write("\n")
    

def main():
    """
    compare SNOTEL station data to WRF Precip
    """

    print("Reading snotel list")
    snotel_stations=load_snotel_list(snotel_list_filename)
    print("Loading WRF Geographic data")
    wrfgeo=load_wrf_geo(wrf_geo_file)
    print("Loading WRF precip data")
    fullwrfdata=ncio.read_nc(wrf_data_file,"prec").data
    stats=[]
    print("looping through stations:")
    
    for station in snotel_stations:
        snoteldata=load_snotel_data(station)
        if type(snoteldata)!=type(-1):
            if (snoteldata.shape[0]>(365*5)):
                wrfdata=load_wrf_data(station,wrfgeo,wrf_data_file,data=fullwrfdata)
                snotelstats=calc_stats(snoteldata[:,-1])
                wrfstats=calc_stats(wrfdata)
                print(station.name)
                stats.append(Bunch(station=station,snotel=snotelstats,wrf=wrfstats))
    print("Writing data")
    write_stats(stats)
    
        
        
if __name__ == '__main__':
    main()