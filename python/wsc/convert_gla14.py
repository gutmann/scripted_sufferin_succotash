#!/usr/bin/env python
import os
import re
import glob
from multiprocessing import Pool

import numpy as np

import swim_io as io

template_file="control.template"
control_file="read_glas_ctrl.dat"
latmin=35.
latmax=45.
lonmin=360-120.
lonmax=360-100.
escale=1000.0
llscale=1000000.0

def write_idl_input(filename):
    idlinputfile=filename.replace(".DAT",".idlinput")
    idlin=open(idlinputfile,"wu")
    idlin.write('!path=+".:"+!path\n')
    idlin.write(".run read_glas_file\n")
    idlin.write("read_glas_file, ControlFile='"+filename.replace(".DAT",".ctrl")+"'\n")
    idlin.write('print, "Finished Reading GLAS"\n')
    idlin.write("exit\n")
    idlin.close()
    

def setup_idl_file(filename):
    template=open(template_file,"ru")
    control=open(filename.replace(".DAT",".ctrl"),"wu")
    for l in template:
        if re.match(".*INPUT_FILE.*",l):
            prefix=l.split("=")[0]
            control.write(prefix+"= "+filename)
        elif re.match(".*OUTPUT_FILE.*",l):
            prefix=l.split("=")[0]
            control.write(prefix+"= "+filename.replace(".DAT",".out"))
        else:
            control.write(l)
    control.close()
    template.close()
    
    write_idl_input(filename)
            

def convert_glas_file(filename):
    setup_idl_file(filename)
    os.system("idl <"+filename.replace(".DAT",".idlinput")+" >>"+filename.replace(".DAT",".idloutput"))




def read_time(filedata):
    found_time=False
    while not found_time:
        line=filedata.next()
        if re.match(".*i_UTCTime.*",line):
            found_time=True
    return np.array([np.float64(line.split()[2])])

def read_lats(filedata):
    found_lats=False
    while not found_lats:
        line=filedata.next()
        if re.match(".*i_lat.*",line):
            found_lats=True
    data=line.split()[2:]
    for i in range(6):
        line=filedata.next()
        data.extend(line.split())
        
    return np.array(data,dtype=np.float64)

def read_lons(filedata):
    line=filedata.next()
    data=line.split()[2:]
    for i in range(6):
        line=filedata.next()
        data.extend(line.split())
        
    return np.array(data,dtype=np.float64)

def read_elev(filedata):
    line=filedata.next()
    data=line.split()[2:]
    for i in range(6):
        line=filedata.next()
        data.extend(line.split())
        
    return np.array(data,dtype=np.float64)


def load_next(filedata):
    bad_data=(np.array([-99]),np.array([-99]),np.array([-99]),np.array([-99]))
    time=read_time(filedata)
    lats=read_lats(filedata)
    if lats.size!=40:
        return bad_data
    lons=read_lons(filedata)
    if lons.size!=40:
        return bad_data
    elevs=read_elev(filedata)
    if elevs.size!=40:
        return bad_data
    time=time.repeat(lats.size)
    return lats/llscale,lons/llscale,time,elevs/escale
    
def read_data(filename):
    lats=[]
    lons=[]
    elevs=[]
    times=[]
    laterr=[]
    lonerr=[]
    fdata=open(filename,"ru")
    moredata=True
    try:
        while moredata:
            curlat,curlon,curtime,curelev=load_next(fdata)
            if curlat[0]==-99:
                print("Missed data:"+fdata.next())
            else:
                tmp=np.where((curlat<latmax)&(curlat>latmin)&(curlon<lonmax)&(curlon>lonmin))
                if len(tmp[0])>0:
                    lats.append(curlat[tmp])
                    lons.append(curlon[tmp])
                    times.append(curtime[tmp])
                    elevs.append(curelev[tmp])
                else:
                    # pass
                    laterr.append([curlat[curlat<1000].min(),curlat[curlat<1000].max()])
                    lonerr.append([curlon[curlon<1000].max(),curlon[curlon<1000].min()])
    except StopIteration:
        fdata.close()
    if len(lats)==0:
        return np.array([-9999]),laterr,lonerr
    lats=np.concatenate(lats)
    lons=np.concatenate(lons)
    times=np.concatenate(times)
    elevs=np.concatenate(elevs)
    
    length=lats.size
    
    outputdata=np.zeros((length,4))
    outputdata[:,0]=lats
    outputdata[:,1]=lons
    outputdata[:,2]=times
    outputdata[:,3]=elevs
    
    return outputdata,laterr,lonerr
    
        

def subset_data(filename):
    outputdata,laterr,lonerr=read_data(filename)
    if outputdata.size>3:
        io.write(filename.replace(".out",".nc"),outputdata,"d")
    
def convert_and_subset_data(filename):
    convert_glas_file(filename)
    subset_data(filename.replace(".DAT",".out"))

def main():
    files=glob.glob("../*.DAT")
    p=Pool(20)
    p.map(convert_and_subset_data,files)
    

if __name__ == '__main__':
    main()

## test
# import glob
# from wsc import convert_gla14
# filename=glob.glob("*.out")[0]
# data=convert_gla14.read_data(filename)
