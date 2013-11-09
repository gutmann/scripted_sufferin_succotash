#!/usr/bin/env python

import gc,sys,glob

import numpy as np
import matplotlib.pyplot as plt

import swim_io

def load_obs():
    return np.concatenate(swim_io.read_files("obs/*","pr"))

def load_loca():
    data=swim_io.read_nc("dp_regridded.pr.historical.1976-2005.yearly.day_cat.srs.bc_presrat_mon.units.ds.nc","pr").data
    return data[:,9:-1,10:-8]

def scaletohuc(data,huc):
    huclist=huc[0]
    hucdex=huc[1]
    outputdata=np.zeros((data.shape[0],len(huclist)))
    for i,h in enumerate(huclist):
        index=np.where((hucdex==h)&(data[0,...]<1e10))
        if len(index[0])>0:
            curdata=data[:,index[0],index[1]]
            # print(h,curdata.shape,len(index[0]),len(index[1]))
            outputdata[:,i]=curdata.mean(axis=1)
            
    return outputdata
    
def load_hucs():
    bd="/d2/gutmann/usbr/hucdata/"
    hucfiles=[bd+"HUC02/huc2_",bd+"HUC04/huc4_",bd+"HUC08/huc8_"]
    hucs=[]
    for h in hucfiles:
        filename=glob.glob(h+"*12km*.nc")[0]
        hucdata=swim_io.read_nc(filename).data
        huclist=np.unique(hucdata)
        hucname=filename.split("/")[-1][:4]
        hucs.append([huclist,hucdata,hucname])
    return hucs
    
    
def calc_wetfrac(data):
    wet=np.zeros(data.shape)
    wet[data>0]=1
    wf=wet.sum(axis=0)/float(wet.shape[0])
    del wet
    gc.collect()
    return wf
    
def calc_p99(data):
    """Calculate the 99th and 1st percentile for data along the first index
    
    data can be a 2D or 3D array, 
        sorts data on axis 0 and 
        returns the 99th and 1st percentile on a 1D or 2D grid
    """
        
    output=np.zeros(data.shape[1:])
    
    for i in range(data.shape[1]):
        print(i, data.shape[1])
        sys.stdout.flush()
        if len(data.shape)>2:
            sorted_data=np.sort(data[:,i,:],axis=0)
            output[i,:]=sorted_data[np.round(data.shape[0]*0.99),...]
        else:
            sorted_data=np.sort(data[:,i])
            output[i]=sorted_data[np.round(data.shape[0]*0.99)]

    return output


def process_data(data,name,hucs):
    print("Processing "+name)
    print("Calculating wetfrac")
    sys.stdout.flush()
    # swim_io.write(name+"_wetfrac_"+"grid",calc_wetfrac(data))
    print("Calculating p99")
    sys.stdout.flush()
    # swim_io.write(name+"_p99_"+"grid",calc_p99(data))
    gc.collect()
    
    for h in hucs:
        print(h[2])
        print("Scaling...")
        sys.stdout.flush()
        hucdata=scaletohuc(data,h)
        print("Calculating wetfrac")
        sys.stdout.flush()
        swim_io.write(name+"_wetfrac_"+h[2],calc_wetfrac(hucdata))
        print("Calculating p99")
        sys.stdout.flush()
        swim_io.write(name+"_p99_"+h[2],calc_p99(hucdata))
    
    
def main():
    hucs=load_hucs()
    
    # print("loading obs")
    # obs=load_obs()
    # process_data(obs,"obs",hucs)
    # del obs
    # gc.collect()
    
    print("loading loca")
    sys.stdout.flush()
    loca=load_loca()
    process_data(loca,"loca",hucs)
    gc.collect()


if __name__ == '__main__':
    main()