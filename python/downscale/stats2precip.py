#!/usr/bin/env python
from __future__ import print_function
import glob,datetime,sys
import numpy as np
from scipy.stats import norm
import swim_io
import nc_reader

output_file="probabilistic_precip"
prcpfile="grid_prcp_norm.nc"
probfile="grid_prob_prcp_norm.nc"
randfilesearch="stochastic/sprand*.nc"
rand_files=glob.glob(randfilesearch)
rand_files.sort()

base_date=datetime.datetime(1970,01,01)
file_variable = ["tmp_2m","dswrf_sfc","cape_sfc","pwat_eatm"]
variable_name = ["TMP_ens_mean_2maboveground",
                 "DSWRF_ens_mean_surface",
                 "CAPE_ens_mean_surface",
                 "PWAT_ens_mean_entireatmosphere_consideredasasinglelayer_"]

def load_rand(i):
    return swim_io.read_nc(rand_files[i]).data

def load_gefs(date,geo_file,geolut=None):
    year=str(date.year)
    month="{0:02}".format(date.month)
    day="{0:02}".format(date.day)
    output=[]
    for f,v in zip(file_variable,variable_name):
        curfile="gefs/"+year+month+day+"/"+f+"*_mean.nc"
        gefs_time=2 #third hour(?) in the file?
        # takes care of all the hard work of geo-interpolating
        # ideally should pass the geolut back out so I don't have to recreate it every time...
        gefs=nc_reader.NC_Reader(curfile,geomatch_file=geo_file,
            glatvar="latitude",glonvar="longitude",readvars=[v],
            ntimes=2,firstfile_timeinit=gefs_time)
        output.append(gefs.next()[0])
        
    return output
    

def main():
    """convert probability coefficients to precip amounts"""
    timeseconds=swim_io.read_nc(probfile,"time").data
    time=[base_date+datetime.timedelta(i/86400.0) for i in timeseconds]
    lat=swim_io.read_nc(probfile,"latitude").data
    lon=swim_io.read_nc(probfile,"longitude").data
    prob=swim_io.read_nc(probfile,"coefficient").data
    prcp=swim_io.read_nc(prcpfile,"coefficient").data
    nx=prob.shape[-1]
    ny=prob.shape[-2]
    nt=prob.shape[-3]
    nv=prob.shape[0]
    nt=40
    
    print(time[0],time[-1])
    output_data=np.zeros((nt,ny,nx))
    output_data2=np.zeros((nt,ny,nx))
    output_data3=np.zeros((nt,ny,nx))
    print("starting")
    for i in range(nt):
        print(i,nt,end="\r")
        sys.stdout.flush()
        gefs_data=load_gefs(time[i],geo_file=probfile)
        rand_data=load_rand(i)
        curprec=prcp[0,i,...]
        curprob=prob[0,i,...]
        for v in range(1,nv):
            curprec+=gefs_data[v-1]*prcp[v,i,...]
            curprob+=gefs_data[v-1]*prob[v,i,...]
        
        # curprec[norm.cdf(rand_data)<curprob]=0
        output_data[i,...]=curprec
        output_data2[i,...]=curprob
        curprob[curprob<norm.cdf(rand_data)]=0
        output_data3[i,...]=curprob
    
    swim_io.write(output_file,output_data)
    swim_io.write(output_file+"prob",output_data2)
    swim_io.write(output_file+"thresh",output_data3)
            
    
    
    

if __name__ == '__main__':
    main()