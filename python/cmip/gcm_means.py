#!/usr/bin/env python
from __future__ import print_function
import sys,os,glob

import numpy as np

import mygis
import units
from cmip import ccsm,gfdl,ipsl,fgoals,cal,access

# used in 3D means
varnames=["hus","va","ua","ta"]

global_filesearch="{varname}_6hrLev_*_historical_r*i*p*_*"

# dictionaries to pull out GCM specific handlers
# defines functions to be used in computing the vertical coordinate.  
# This computes pressure for all but access which computes geopotential height and pressure has to be infered...?
global_vert_coords=dict(
                        access=access.vcoord,
                        bcc = ccsm.vcoord,
                        bnu = ccsm.vcoord,
                        canesm = ipsl.vcoord,
                        ccsm = ccsm.vcoord,
                        cnrm_cm5 = ccsm.vcoord,
                        fgoals = fgoals.vcoord,
                        gfdl_cm3 = gfdl.vcoord,
                        gfdl_esm = gfdl.vcoord,
                        ipsl_mr = ipsl.vcoord,
                        miroc5 = ccsm.vcoord,
                        miroc_esm = ccsm.vcoord,
                        mk3 = ccsm.vcoord,
                        mri_cgcm3 = ccsm.vcoord,
                        noresm = ccsm.vcoord,
                        giss_e2h = ccsm.vcoord
                        )
# defines the functions to use for the calendar system for each model.  
# cal is a calendar module, std_data = gregorian
global_timeing=dict(
                    access=cal.std_date,
                    bcc = cal.noleap_date,
                    bnu = cal.noleap_date,
                    canesm = cal.noleap_date,
                    ccsm=cal.noleap_date,
                    cnrm_cm5 = cal.std_date,
                    fgoals = cal.noleap_date,
                    gfdl_cm3 = cal.noleap_date,
                    gfdl_esm = cal.noleap_date,
                    ipsl_mr = cal.noleap_date,
                    miroc5 = cal.noleap_date,
                    miroc_esm = cal.std_date,
                    mk3 = cal.noleap_date,
                    mri_cgcm3 = cal.std_date,
                    noresm = cal.noleap_date,
                    giss_e2h = cal.noleap_date,
                    )

# this is the year which the calendar system for each GCM is based off of.  
# None means it changes from file to file and needs to be read explicitly
# all could be set to None and it would probably work just fine...
global_y0=dict(
                bcc = None,
                bnu = 1950,
                canesm = 1850,
                ccsm=1850,
                cnrm_cm5 = 1850,
                fgoals = 1,
                gfdl_cm3 = 1860,
                gfdl_esm = 1861,
                ipsl_mr = 1850,
                miroc5 = 1850,
                miroc_esm = 1850,
                mk3 = 1850,
                mri_cgcm3 = 1850,
                noresm = 1850,
                giss_e2h = None,
                access = 1
                )
# like y0 but for starting month again all could be set to None and it should work
global_m0=dict(
                bcc =       1,
                bnu =       1,
                canesm =    1,
                ccsm=       1,
                cnrm_cm5 =  1,
                fgoals =    1,
                gfdl_cm3 =  1,
                gfdl_esm =  1,
                ipsl_mr =   1,
                miroc5 =    1,
                miroc_esm = 1,
                mk3 =       1,
                mri_cgcm3 = 1,
                noresm =    1,
                giss_e2h =    None,
                access =    1
                )

# The years to compute climatological means over (inclusive)
#  start in 1979 because that is the first year for most high quality reanalysis data. 
start_year=1979
end_year=2005

def med_filt(data,minval=0):
    outputdata=np.zeros(data.shape)
    for i in range(data.shape[0]-2):
        for j in range(data.shape[1]-2):
            window=data[i:i+3,j:j+3]
            if len(np.where(window>minval)[0])>0:
                outputdata[i+1,j+1]=np.median(window[window>minval])
    outputdata[0,:] =data[0,:]
    outputdata[-1,:]=data[-1,:]
    outputdata[:,0] =data[:,0]
    outputdata[:,-1] =data[:,-1]
    return outputdata



def calc_mean_pressure_levels(gcm="ccsm"):
    calc_date=global_timeing[gcm]
    y0=global_y0[gcm]
    m0=global_m0[gcm]
    calc_vert_coord=global_vert_coords[gcm]
    filesearch=global_filesearch
    
    output_data=[dict() for i in range(12)]
    
    files=glob.glob(filesearch.format(varname=varnames[0]))
    files.sort()
    for f in files:
        # start with an arbitrary (must be 3D) variable to read timing and 3D coordinates
        base_file=f
        dates=mygis.read_nc(base_file,"time").data
        last_date=calc_date(dates[-1],y0=y0,m0=m0,filename=base_file)
        print(last_date.year)
        # print("{0.year}/{0.month}/{0.day} {0.hour}:{0.minute}:{0.second}".format(last_date))
        if last_date.year>=start_year:
            pressures=calc_vert_coord(base_file)
            for i,d in enumerate(dates):
                curdate=calc_date(d,y0=y0,m0=m0,filename=base_file)
                if curdate.year>=start_year:
                    if "p" in output_data[curdate.month-1]:
                        output_data[curdate.month-1]["p"]+=pressures[i,...]
                        output_data[curdate.month-1]["p_n"]+=1
                    else:
                        output_data[curdate.month-1]["p"]=np.zeros(pressures[i,...].shape)
                        output_data[curdate.month-1]["p"][:]=pressures[i,...]
                        output_data[curdate.month-1]["p_n"]=1
    for i in range(12):
        output_data[i]["p"]/=output_data[i]["p_n"]
        mygis.write(gcm+"_month{0:02}_mean_plevel.nc".format(i+1),data=output_data[i]["p"],varname="pressure")

def calc_2d_means(gcm,varname="ps",filesearch=None):
    """docstring for calc_2d_means"""
    calc_date=global_timeing[gcm]
    y0=global_y0[gcm]
    if varname=="sst":
        y0=None
    m0=global_m0[gcm]
    if (filesearch==None):
        filesearch=global_filesearch.format(varname=varname)
    
    output_data=[dict() for i in range(12)]
    
    files=glob.glob(filesearch)
    files.sort()
    for f in files:
        print(f)
        # start with an arbitrary (must be 3D) variable to read timing and 3D coordinates
        base_file=f
        dates=mygis.read_nc(base_file,"time").data
        last_date=calc_date(dates[-1],y0=y0,m0=m0,filename=base_file)
        if last_date.year>=start_year:
            print("{year}/{month}/{day} {hour}:{minute}:{second}".format(**last_date))
            
            curdata=mygis.read_nc(f,varname).data
            for i,d in enumerate(dates):
                curdate=calc_date(d,y0=y0,m0=m0,filename=base_file)
                if curdate.year>=start_year:
                    if varname=="sst":
                        #median filter pushes data to cover more land, and removes random speckling in some datasets (e.g. CNRM!)
                        curdata[i,...]=med_filt(curdata[i,...])
                    
                    if varname in output_data[curdate.month-1]:
                        output_data[curdate.month-1][varname]+=curdata[i,...]
                        output_data[curdate.month-1][varname+"_n"]+=1
                    else:
                        output_data[curdate.month-1][varname]=np.zeros(curdata[i,...].shape)
                        output_data[curdate.month-1][varname][:]=curdata[i,...]
                        output_data[curdate.month-1][varname+"_n"]=1
    for i in range(12):
        output_data[i][varname]/=output_data[i][varname+"_n"]
        mygis.write(gcm+"_month{0:02}_mean_{1}.nc".format(i+1,varname),data=output_data[i][varname],varname=varname)

def load_orography(gcm):
    landorog = mygis.read_nc("orography.nc","orog").data
    return landorog
    # landfrac = mygis.read_nc("sftlf.nc","sftlf").data
    # return landorog*(landfrac/100.0)

def compute_3d_means(gcm):
    """docstring for compute_all_means"""
    calc_date=global_timeing[gcm]
    y0=global_y0[gcm]
    m0=global_m0[gcm]
    calc_vert_coord=global_vert_coords[gcm]
    try:
        z=load_orography(gcm)
    except:
        print("Could not load orography")
        z=None
    
    filesearch=global_filesearch
    varfiles=dict()
    for v in varnames:
        curfiles=glob.glob(filesearch.format(varname=v))
        curfiles.sort()
        varfiles[v]=curfiles
    nfiles=len(curfiles)
    
    output_data=[dict() for i in range(12)]
    
    for filenum in range(nfiles):
        base_file=varfiles[varnames[0]][filenum]
        dates=mygis.read_nc(base_file,"time").data
        
        print(calc_date(dates[-1],y0=y0,m0=m0,filename=base_file).year,end="\r")
        sys.stdout.flush()
        if (calc_date(dates[-1],y0=y0,m0=m0,filename=base_file).year>=start_year):
            cur_pressures=calc_vert_coord(base_file)
            curdata=dict(hus=None,ta=None,slp=None,psl=None,ps=None)
            for v in varnames:
                if len(varfiles[v])>filenum:
                    print(varfiles[v][filenum])
                    curdata[v]=mygis.read_nc(varfiles[v][filenum],v).data
                else:
                    print(len(varfiles[v]),filenum,varfiles[v][-1])
                    curdata[v]=None
                    
            if np.max(cur_pressures)<90000:
                print(np.max(cur_pressures))
                # cur pressures is actually geopotential height (in ACCESS 1.3 GCM)
                # so we need to calculated the 3D pressure field
                print("Computing Pressures : WARNING this was designed for the ACCESS GCM, it may work for others, but needs to be tested")
                print("   Assumes 10 3D atm files per surface pressure file among other things")
                #cur pressures is actually geopotential height
                if curdata["ps"]==None:
                    psfiles=glob.glob(filesearch.format(varname="ps"))
                    psfiles.sort()
                    print(psfiles[filenum/10])
                    ps=mygis.read_nc(psfiles[filenum/10],"ps").data
                    pstime=mygis.read_nc(psfiles[filenum/10],"time").data
                    
                ps_timestep=np.where(pstime==dates[0])[0]
                levels=cur_pressures.repeat(dates.size,axis=0)
                curdata["ps"]=ps[ps_timestep:ps_timestep+curdata["ta"].shape[0]]
                cur_pressures=units.zt2p(levels-z[np.newaxis,np.newaxis,...],
                                         p0=curdata["ps"],
                                         t0=curdata["ta"])
                cur_pressures[~np.isfinite(cur_pressures)]=0.1
                cur_pressures[cur_pressures<0.1]=0.1
            else:
                levels=None
            print("Computing RH")
            if (curdata["ta"]!=None) and (curdata["hus"]!=None):
                curdata["rh"]=units.sh2rh(curdata["ta"],cur_pressures/100.0,curdata["hus"])
            else:
                curdata["rh"]=None
            
            
            print("Computing SLP")
            if (curdata["ta"]!=None) and (curdata["hus"]!=None) and (z!=None):
                try:
                    if gcm=="access":
                        ps=curdata["ps"]
                    else:
                        ps=mygis.read_nc(base_file,"ps").data
                    
                    curdata["slp"]=units.calc_slp(ps,z[np.newaxis,...],ts=curdata["ta"][:,0,...],mr=curdata["hus"][:,0,...],method=2)
                except Exception as e:
                    print(e)
                    print("Could not compute Sea level pressure. ")
            else:
                if curdata["psl"]!=None:
                    curdata["slp"]=curdata["psl"]

            
            if (levels!=None):
                curdata["z"]=levels
            else:
                print("Computing Geopotential Height")
                if (curdata["ta"]!=None) and (curdata["hus"]!=None) and (curdata["slp"]!=None):
                    try:
                        curdata["z"]=np.zeros(curdata["ta"].shape)
                        for z_time in range(curdata["ta"].shape[0]):
                            curdata["z"][z_time,...]=units.calc_z(curdata["slp"][z_time],
                                                                  cur_pressures[z_time],
                                                                  t=curdata["ta"][z_time],
                                                                  mr=curdata["hus"][z_time])
                    except Exception as e:
                        print(e)
                        print("Could not compute Geopotential height. ")
            
            print("Storing data...")
            for i in range(len(dates)):
                curdate=calc_date(dates[i],y0=y0,m0=m0,filename=base_file)
                if (curdate.year>=start_year) and (curdate.year<=end_year):
                    for v in curdata.keys():
                        if curdata[v]!=None:
                            if v in output_data[curdate.month-1]:
                                output_data[curdate.month-1][v]+=curdata[v][i,...]
                                output_data[curdate.month-1][v+"_n"]+=1
                            else:
                                output_data[curdate.month-1][v]=np.zeros(curdata[v][i,...].shape,dtype=np.float64)
                                output_data[curdate.month-1][v][:]=curdata[v][i,...]
                                output_data[curdate.month-1][v+"_n"]=1
    
    print("Writing data to files")
    for i in range(12):
        for v in curdata.keys():
            print(i,v)
            if v in output_data[i].keys():
                output_filename=gcm+"_month{0:02}_mean_{1}.nc".format(i+1,v)
                print(output_filename)
                mygis.write(output_filename,
                            data=output_data[i][v]/output_data[i][v+"_n"],varname=v)
            else:
                print("ERROR: variable: {} for month {} does not exist".format(v,str(i)))
    
    
def process_all_gcms():
    """docstring for process_all_gcms"""
    cwd=os.getcwd()
    for gcm in global_vert_coords.keys():
        print(gcm)
        os.chdir(cwd)
        try:
            os.chdir(gcm+"/subset")
            
            calc_2d_means(gcm,"ps")
            calc_2d_means(gcm,"psl",filesearch="psl_6hrPlev_*_historical_r*i*p*_*.nc")
            # calc_mean_pressure_levels(gcm)
            
            compute_3d_means(gcm)
            calc_2d_means(gcm,"sst",filesearch="sst_day_*historical_r*i*p*_*.nc")
            
        except Exception as e:
            print("Failed to process: "+gcm)
            print(e)

def main(gcm="ccsm",pressure_only=True):
    """docstring for main"""
    print("Computing historical monthly means for model: "+gcm)
    if not pressure_only:
        # pressure_levels=mygis.read_files(gcm+"_month*_plevel.nc","pressure")
        compute_3d_means(gcm)
        calc_2d_means(gcm,"ps")
        calc_2d_means(gcm,"psl",filesearch="psl_6hrPlev_*_historical_r*i*p*_*.nc")
        calc_2d_means(gcm,"sst",filesearch="sst_day_*historical_r*i*p*_*.nc")
        
    else:
        calc_mean_pressure_levels(gcm)
    

if __name__ == '__main__':
    pressure_only=False
    gcm="ccsm"
    if len(sys.argv)>1:
        gcm=sys.argv[1]
    if len(sys.argv)>2:
        pressure_only=(sys.argv[2]=="pressure_only")
    
    if gcm=="batch":
        process_all_gcms()
    else:
        main(gcm=gcm,pressure_only=pressure_only)
    # try:
    #     main(gcm=gcm,pressure_only=pressure_only)
    # except KeyError as e:
    #     print("Key Error: "+str(e))
    #     print("May not know how to run for model:"+gcm+" yet.")
    