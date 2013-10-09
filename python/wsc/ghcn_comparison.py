#!/usr/bin/env python
import glob,re
import numpy as np

from bunch import Bunch
import swim_io as ncio


ghcn_list_filename="ghcn/ghcn_site_list.txt"
wrf_geo_file="wrf/4km_wrf_output.nc"
wrf_data_file="wrf/wrf_prcp_TIMEERROR.nc"

def load_ghcn_list(ghcn_list_file):
    """Load station name/location information from a text file
    
    example file:
    NSITES 1208
    #STNID            LAT             LONG      ELV_M     STN_NAME
    US1COAD0001,    39.9460000,    -104.6088000,  1576, "US1COAD0001"
    USC00020482,    35.3000000,    -112.4833000,  1617, "USC00020482"
    """
    with open(ghcn_list_file) as f:
        for i in range(2):l=f.next()
        
        outputlist=[]
        for l in f:
            line=l.split(",")
            name=line[0].strip()
            lat=np.double(line[1])
            lon=np.double(line[2])
            outputlist.append(Bunch(name=name,lat=lat,lon=lon))
        
    return outputlist

def load_ghcn_data(ghcn_station):
    """load data from a GHCN file 
    file format from andy newman (originally Kyoko?)
    
    #DATE/TIME    = YYYYMMDD HHMMSS
    #PRCP         = Total Precipitation last 24 hours (mm)
    # Missing data     = -999
    #

      DATE   HHMMSS  PRCP
    19800101 080000     -999.0
    19800102 080000     -999.0
    """
    try:
        filename=glob.glob("ghcn/"+ghcn_station.name+"*PRCP*.txt")[0]
        with open(filename) as f:
            l=""
            while (not re.match(".*DATE .*HHMMSS .*PRCP.*",l)):l=f.next()
        
            outputlist=[]
            for l in f:
                line=l.split()
                year=int(l[:4])
                month=int(l[4:6])
                day=int(l[6:8])
            
                precip=float(line[2])
                outputlist.append([year,month,day,precip])
        return np.array(outputlist)
        
    except Exception as e:
        print("ERROR:"+e)
        print("  With station:"+ghcn_station.name)
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
    with open("output_statistics2.txt","w") as f:
        f.write("lat lon ")
        for k in stats[0].wrf.keys():
            f.write(k+" ")
        f.write("\n")
        for stat in stats:
            f.write(str(stat.station.lat)+" "+str(stat.station.lon))
            for k in stat.wrf.keys():
                f.write(" ")
                f.write(str(stat.ghcn[k])+" "+str(stat.wrf[k]))
            f.write("\n")
    

def main():
    """
    compare GHCN station data to WRF Precip
    """

    print("Reading GHCN list")
    ghcn_stations=load_ghcn_list(ghcn_list_filename)
    print("Loading WRF Geographic data")
    wrfgeo=load_wrf_geo(wrf_geo_file)
    print("Loading WRF precip data")
    fullwrfdata=ncio.read_nc(wrf_data_file,"prec").data
    stats=[]
    print("looping through stations:")
    
    for station in ghcn_stations:
        ghcndata=load_ghcn_data(station)
        if type(ghcndata)!=type(-1):
            if (ghcndata.shape[0]>(365*5)):
                wrfdata=load_wrf_data(station,wrfgeo,wrf_data_file,data=fullwrfdata)
                ghcnstats=calc_stats(ghcndata[:,-1])
                wrfstats=calc_stats(wrfdata)
                print(station.name)
                stats.append(Bunch(station=station,ghcn=ghcnstats,wrf=wrfstats))
    print("Writing data")
    write_stats(stats)
    
        
        
if __name__ == '__main__':
    main()