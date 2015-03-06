#!/usr/bin/env python
import glob
import sys

import numpy as np
import mygis

outputfile="means_{}.txt"
precip_var="Prec"
tmax_var="Tmax"
tmin_var="Tmin"

def usage():
    """docstring for usage"""
    print("""USAGE:
    extract_timeseries.py <lat> <lon> <name> [variable]
    
        lat = latitude location to extract
        lon = longitude location to extract
        name= name of station (for output filename)
        variable = name of the variable in *.nc files to extract 
            (defaults to "Prec")
    """)

def find_point(lat_point,lon_point,filename):
    """docstring for calc_point"""
    geo=mygis.read_geo(filename)
    dists=((geo.lat-lat_point)**2) + ((geo.lon-lon_point)**2)
    
    point=np.unravel_index(np.argmin(dists),dists.shape)
    print(lat_point,geo.lat[point[0],point[1]])
    print(lon_point,geo.lon[point[0],point[1]])
    return point

def load_point(files,point):
    """docstring for load_point"""
    outputdata=[]
    year=int(files[0].split(".")[-2][:4])
    year+=int(files[0].split(".")[-2][4]) # convert to water year
    precip=0
    temperature=0
    n=0
    
    for f in files:
        curyear=int(f.split(".")[-2][:4])
        curyear+=int(f.split(".")[-2][4]) # convert to water year
        if year!=curyear:
            temperature/=float(n)
            outputdata.append([year,precip,temperature])
            year=curyear
            precip=0
            temperature=0
            n=0
        
        # load the precipitation data (summed over the year)
        ncdata=mygis.read_nc(f,precip_var,returnNCvar=True)
        tmpp=np.sum(ncdata.data[:,point[0],point[1]])
        precip+=tmpp
        ncdata.ncfile.close()
        
        # load the tmax data (summed over the year then divided by n)
        ncdata=mygis.read_nc(f,tmax_var,returnNCvar=True)
        tmp=ncdata.data[:,point[0],point[1]]
        temperature+=np.sum(tmp)
        n+=ncdata.data.shape[0]
        ncdata.ncfile.close()
        
        # load the tmin data (summed over the year then divided by n)
        ncdata=mygis.read_nc(f,tmin_var,returnNCvar=True)
        tmp2=ncdata.data[:,point[0],point[1]]
        temperature+=np.sum(tmp2)
        n+=ncdata.data.shape[0]
        ncdata.ncfile.close()
        
        print(f.split(".")[-2]+"  "+str(tmpp)+"  "+str(np.mean(tmp)/2+np.mean(tmp2)/2))
    
    
    if n>350:
        temperature/=float(n)
        outputdata.append([year,precip,temperature])
    return np.array(outputdata)

def write_data(data,filename):
    """docstring for write_data"""
    np.savetxt(filename,data,
                fmt=['%i',"%f","%f"])
                # header="Year  Precip[mm]  Temperature[C]",
                
def main(args=[]):
    """docstring for main"""
    if len(args)==0:
        args=sys.argv[1:]
    if len(args)<3:
        usage()
        sys.exit()
    
    lat=float(args[0])
    lon=float(args[1])
    name=args[2]
    
    files=glob.glob("*.nc")
    files.sort()
    xy=find_point(lat,lon,files[0])
    print("location: lat:{}  lon:{}".format(lat,lon))
    print("grid pnt:   y:{}    x:{}".format(xy[0],xy[1]))
    
    data=load_point(files,point=xy)
    
    write_data(data,outputfile.format(name))

if __name__ == '__main__':
    main()