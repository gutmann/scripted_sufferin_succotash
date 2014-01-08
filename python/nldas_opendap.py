#!/usr/bin/env python
# from __future__ import print_function
import numpy as np
from netCDF4 import Dataset
import date_fun

dataset="LSM"
#dataset="Forcing"

def load_coords(filename="station_coords"):
    outputlist=[]
    with open(filename,"ru") as f:
        for l in f:
            curline=l.split()
            name=curline[0]
            lat=np.double(curline[1])
            lon=np.double(curline[2])
            if (lat<=53) and (lat>=25) and (lon<=-67) and (lon>=-125):
                outputlist.append([name,lat,lon])
            else:
                print(name+" outside NLDAS Bounds")
    return outputlist
            

def main(site_lat=41.0674,site_lon=-106.1218,outputfile="stationID.txt",
         varlist=["ugrd10m","vgrd10m"],
         nldas_url="http://hydro1.sci.gsfc.nasa.gov/dods/NLDAS_FORA0125_H.002"):
        
    nldas_to=date_fun.datearr2mjd(np.array([1,1,1,0,0,0.0]))
    begintime=date_fun.datearr2mjd(np.array([2013,1,1,0,0,0.0]))-nldas_to
    endtime=date_fun.datearr2mjd(np.array([2014,1,1,0,0,0.0]))-nldas_to

    nldas=Dataset(nldas_url)
    times=nldas.variables["time"][:]
    lat=nldas.variables["lat"][:]
    lon=nldas.variables["lon"][:]
    firsttime=np.where(times>=begintime)[0][0]
    lasttime=np.where(times>=endtime)[0][0]
    
    x=np.argmin(np.abs(lon-site_lon))
    y=np.argmin(np.abs(lat-site_lat))
        
    outputdata=[]
    for v in varlist:
        print("   "+v)
        outputdata.append(nldas.variables[v][firsttime:lasttime,y,x])

    for (i,time) in enumerate(range(firsttime,lasttime)):
        date=date_fun.mjd2date(times[time]+nldas_to,roundseconds=True)
    
    fout=open(outputfile,"w")
    fout.write("Year Month Day Hour Minute ")
    for v in varlist:
        fout.write("   "+v)
    fout.write("\n")
    
    for i in range(lasttime-firsttime):
        date=date_fun.mjd2date(times[i+firsttime]+nldas_to,roundseconds=True)
        fout.write(str(date[:5]).replace('[','').replace(']','  '))
        for v in outputdata:
            fout.write(" "+str(v[i]))
        fout.write("\n")
    fout.close()
    
    

if __name__ == '__main__':
    if dataset=="LSM":
        lsm_url="http://hydro1.sci.gsfc.nasa.gov/dods/NLDAS_NOAH0125_H.002"
        varlist=["soilm0_10cm","soilm10_40cm","soilm40_100cm","soilm100_200cm","tsoil0_10cm","lhtflsfc","shtflsfc","snodsfc"]
        # varlist=["soilm0_10cm","tsoil0_10cm"]
       # main(varlist=varlist,nldas_url=lsm_url)
                
        sitelist=load_coords()
        for site in sitelist:
            print(site[0]+","+str(site[1])+","+str(site[2]))
            main(site_lat=site[1],site_lon=site[2],outputfile=site[0]+".txt",
                   varlist=varlist,nldas_url=lsm_url)
            print("...Completed\n")
    else:
        main()
