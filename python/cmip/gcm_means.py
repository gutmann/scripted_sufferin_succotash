#!/usr/bin/env python
import sys
import glob

import numpy as np

import mygis
import units
from cmip import ccsm,gfdl,ipsl,fgoals,cal

varnames=["hus","va","ua","ta"]
translate=dict(hus="RH",va="va",ua="ua",ta="ta",ps="ps")

# dictionaries to pull out GCM specific handlers
global_filesearch="{varname}_6hrLev_*_historical_r?i?p?_*"
global_vert_coords=dict(
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
                        noresm = ccsm.vcoord
                        )
global_timeing=dict(
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
                        noresm = cal.noleap_date
                        )

start_year=1979


def calc_mean_pressure_levels(gcm="ccsm"):
    calc_date=global_timeing[gcm]
    calc_vert_coord=global_vert_coords[gcm]
    filesearch=global_filesearch
    
    output_data=[dict() for i in range(12)]
    
    files=glob.glob(filesearch.format(varname=varnames[0]))
    files.sort()
    for f in files:
        print(f)
        # start with an arbitrary (must be 3D) variable to read timing and 3D coordinates
        base_file=f
        dates=mygis.read_nc(base_file,"time").data
        last_date=calc_date(dates[-1])
        # print("{0.year}/{0.month}/{0.day} {0.hour}:{0.minute}:{0.second}".format(last_date))
        if last_date.year>=start_year:
            pressures=calc_vert_coord(base_file)
            for i,d in enumerate(dates):
                curdate=calc_date(d)
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

def calc_ps_means(gcm):
    """docstring for calc_ps_means"""
    calc_date=global_timeing[gcm]
    calc_vert_coord=global_vert_coords[gcm]
    filesearch=global_filesearch.format(varname="ps")
    
    output_data=[dict() for i in range(12)]
    
    files=glob.glob(filesearch)
    files.sort()
    for f in files:
        print(f)
        # start with an arbitrary (must be 3D) variable to read timing and 3D coordinates
        base_file=f
        dates=mygis.read_nc(base_file,"time").data
        last_date=calc_date(dates[-1])
        # print("{0.year}/{0.month}/{0.day} {0.hour}:{0.minute}:{0.second}".format(last_date))
        if last_date.year>=start_year:
            pressures=mygis.read_nc(f,"ps").data
            for i,d in enumerate(dates):
                curdate=calc_date(d)
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
        mygis.write(gcm+"_month{0:02}_mean_ps.nc".format(i+1),data=output_data[i]["p"],varname="ps")

def compute_3d_means(gcm,pressure_levels):
    """docstring for compute_all_means"""
    calc_date=global_timeing[gcm]
    calc_vert_coord=global_vert_coords[gcm]
    filesearch=global_filesearch
    varfiles=dict()
    for v in varnames:
        curfiles=glob.glob(filesearch.format(varname=v))
        curfiles.sort()
        varfiles[v]=curfiles
    nfiles=len(curfiles)
    output_data=[dict() for i in range(12)]
    
    for filenum in range(nfiles):
        basefile=varfiles[varnames[0]][filenum]
        dates=mygis.read_nc(basefile,"time").data
        cur_pressures=calc_vert_coord(basefile)
        
        if (calc_date(dates[-1]).year>=start_year):
            curdata=dict()
            for v in varnames:
                print(varfiles[v][filenum])
                curdata[v]=mygis.read_nc(varfiles[v][filenum],v).data
            print("Computing RH")
            curdata["rh"]=units.sh2rh(curdata["ta"],cur_pressures/100.0,curdata["hus"])
            print("Storing data...")
            for i in range(len(dates)):
                curdate=calc_date(dates[i])
                if (curdate.year>=start_year):
                    for v in curdata.keys():
                        if v in output_data[curdate.month-1]:
                            output_data[curdate.month-1][v]+=curdata[v][i,...]
                            output_data[curdate.month-1][v+"_n"]+=1
                        else:
                            output_data[curdate.month-1][v]=np.zeros(curdata[v][i,...].shape,dtype=np.float64)
                            output_data[curdate.month-1][v][:]=curdata[v][i,...]
                            output_data[curdate.month-1][v+"_n"]=1
                    
    for i in range(12):
        for v in curdata.keys():
            print(i,v)
            if v in output_data[i].keys():
                output_filename=gcm+"_month{0:02}_mean_{1}.nc".format(i+1,v)
                print(output_filename)
                mygis.write(output_filename,
                            data=output_data[i][v]/output_data[i][v+"_n"],varname=v)
    
    

def main(gcm="ccsm",pressure_only=True):
    """docstring for main"""
    print("Computing historical monthly means for model: "+gcm)
    if not pressure_only:
        # pressure_levels=mygis.read_files(gcm+"_month*_plevel.nc","pressure")
        pressure_levels=None
        compute_3d_means(gcm,pressure_levels)
        # calc_ps_means(gcm)
    else:
        calc_mean_pressure_levels(gcm)
    

if __name__ == '__main__':
    pressure_only=False
    gcm="ccsm"
    if len(sys.argv)>1:
        gcm=sys.argv[1]
    if len(sys.argv)>2:
        pressure_only=(sys.argv[2]=="pressure_only")
        
    try:
        main(gcm=gcm,pressure_only=pressure_only)
    except KeyError as e:
        print("Key Error: "+str(e))
        print("May not know how to run for model:"+gcm+" yet.")
    