#!/usr/bin/env python
from __future__ import print_function
import glob,datetime,sys
import numpy as np
from scipy.stats import norm
import swim_io
import nc_reader

output_file="probabilistic_precip"
prcpfile="grid_04_prcp_norm.nc"
probfile="grid_04_prob_prcp_norm.nc"
randfilesearch="stochastic/sprand*.nc"
rand_files=glob.glob(randfilesearch)
rand_files.sort()

base_date=datetime.datetime(1970,01,01)
nvars=6

if nvars==3:
    file_variable = ["tmp_2m","cape_sfc","pwat_eatm"]
if nvars==4:
    file_variable = ["tmp_2m","cape_sfc","dswrf_sfc","pwat_eatm"]
if nvars==6:
    file_variable = ["tmp_2m","cape_sfc","dswrf_sfc","pwat_eatm","ugrd_80","vgrd_80"]

variable_name = ["TMP_ens_mean_2maboveground",
                 "CAPE_ens_mean_surface",
                 "PWAT_ens_mean_entireatmosphere_consideredasasinglelayer_"]
if nvars>=4:
    variable_name.insert(2,"DSWRF_ens_mean_surface")
if nvars==6:
    variable_name.append("UGRD_ens_mean_80maboveground")
    variable_name.append("VGRD_ens_mean_80maboveground")
    

def load_rand(i):
    return swim_io.read_nc(rand_files[i]).data

def load_gefs(date,geo_file,geolut=None):
    year=str(date.year)
    month="{0:02}".format(date.month)
    day="{0:02}".format(date.day)
    output=[]
    for f,v in zip(file_variable,variable_name):
        curfile="gefs/"+year+month+day+"/"+f+"*_mean.nc"
        gefs_time=7 #21st hour(?) in the file?
        # takes care of all the hard work of geo-interpolating
        # ideally should pass the geolut back out so I don't have to recreate it every time...
        gefs=nc_reader.NC_Reader(curfile,geomatch_file=geo_file,
            glatvar="latitude",glonvar="longitude",readvars=[v],
            ntimes=2,firstfile_timeinit=gefs_time,nn=False,bilin=True)
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
    nt=len(time)
    print(nv)
    print(file_variable,variable_name)
    # nt=30
    print(time[0],time[-1],len(time),nt)
    output_data=np.zeros((nt,ny,nx))
    output_data2=np.zeros((nt,ny,nx))
    output_data3=np.zeros((nt,ny,nx))
    output_data4=np.zeros((nt,ny,nx))
    print("starting")
    for i in range(nt):
        print("  ",i," / ",nt,end="\r")
        sys.stdout.flush()
        gefs_data=load_gefs(time[i],geo_file=probfile)
        rand_data=load_rand(i)
        curprec=prcp[0,i,...]
        curprob=prob[0,i,...]
        for v in range(1,nv):
            # gefs_data[v-1]=gefs_data[v-1]*gains[v-1]+offsets[v-1]
            curprec+=gefs_data[v-1]*prcp[v,i,...]
            curprob+=gefs_data[v-1]*prob[v,i,...]
        
        # curprec[norm.cdf(rand_data)<curprob]=0
        output_data[i,...]=curprec
        output_data2[i,...]=curprob
        curprob[curprob<norm.cdf(rand_data)]=0
        output_data3[i,...]=curprob
        curprec[curprob<norm.cdf(rand_data)]=0
        output_data4[i,...]=curprec**3
    
    swim_io.write(output_file,output_data)
    swim_io.write(output_file+"prob",output_data2)
    swim_io.write(output_file+"prob_thresh",output_data3)
    swim_io.write(output_file+"prec_thresh",output_data4)
            
    
    
    

if __name__ == '__main__':
    main()