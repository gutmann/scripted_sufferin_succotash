#!/usr/bin/env python
import datetime

import matplotlib.pyplot as plt
import numpy as np
import process
import mygis as swim_io

# Smoothing window=870, (=290 when decimation=3000)
# r=0.9918143173048467,
# Polyfit output: array([ -6.09249686e-09,   4.67052267e-04,  -1.10326837e+01, 1.08610472e+05])
# decimation_factor=1000
# start_point=4400*10000
# smooth_min=800
# smooth_max=1000
# smooth_step=10

# 15th order polynomial based off look up table at http://www.adafruit.com/datasheets/103_3950_lookuptable.pdf
# error in fit is <0.1C
temperature_poly=np.array([ -1.84051175e-29,   3.89058030e-26,  -3.72543698e-23,
                         2.13799701e-20,  -8.19978433e-18,   2.21889349e-15,
                        -4.36059768e-13,   6.31196462e-11,  -6.75640933e-09,
                         5.32597237e-07,  -3.05643502e-05,   1.25318268e-03,
                        -3.57800882e-02,   6.92008920e-01,  -9.20463379e+00,
                         7.36758152e+01])

decimation_factor=3000 # convert 10Hz to 5min data)
# start_point=4400*10000l # using full data set
# end_point=None
start_point=1200000 # for a subset of data
master_end_point=0.7 #3600000
# start_point=11000*10000l #!Feb 1
master_smooth_min=1 #200 for freezing?
master_smooth_max=400
master_smooth_step=10

tree=3
if tree==1:
    pot1=3
    pot2=4
if tree==2:
    pot1=5
    pot2=6
if tree==3:
    pot1=12
    pot2=13

def decimate(data,factor):
    """docstring for decimate"""
    n=data.size
    usable=np.floor(float(n)/factor)*factor
    
    data=data[:usable].reshape((usable/factor,factor))
    return data.mean(axis=1)

def smooth_data(data,factor):
    """docstring for smooth_data"""
    newdata=np.zeros(data.shape)
    for i in range(factor,len(data)):
        newdata[i]=data[i-factor:i].mean()
    return newdata

def correct_temperature(data=None,v=None,vo=None,vv=None,r=None,tcol=8,vcol=-1):
    """Calculate the temperature of a thermister from a hard coded polynomial function
    
    Input can either be 
        data: the raw data array from a datalogger with temperature in [:,tcol] and excitation voltage in [:,vcol]
        v and vo (units must match)
            v=thermister voltage (e.g. either voltage or dn)
            vo=excitation voltage (e.g. either voltage of dn)
        vv: v/vo
        or
        r: resistance of the thermister (hard coded as 10*vv/1-vv)
    """
    if r==None:
        if vv==None:
            if vo==None:
                if data==None:
                    raise("ERROR")
                else:
                    vo=data[:,vcol]
                    v=data[:,tcol]
            
            vv=v/vo.astype("d")
        
        r=10*vv/(1-vv)
    temperature=np.zeros(len(r))
    for i in range(len(temperature_poly)):
        temperature+=temperature_poly[len(temperature_poly)-i-1]*r**i
    
    return temperature
    

def optimal_t_fit(data=None,inputT=None,inputCompression=None,dates=None,smin=None,smax=None,sstep=None):
    
    if data==None:
        data=swim_io.read_nc("full_interception.nc").data
    if smin!=None:
        smooth_min=smin
    else:
        smooth_min=master_smooth_min
    if smax!=None:
        smooth_max=smax
    else:
        smooth_max=master_smooth_max
    if sstep!=None:
        smooth_step=sstep
    else:
        smooth_step=master_smooth_step
    
    if master_end_point==None:
        end_point=data.shape[0]
    elif type(master_end_point)==float:
        if master_end_point<1:
            end_point=int(data.shape[0]*master_end_point)
        else:
            end_point=master_end_point
    else:
        end_point=master_end_point
    
    if inputT==None:
        t=decimate(data[start_point:end_point,8].astype("d"),decimation_factor)
        vo=decimate(data[start_point:end_point,-1].astype("d"),decimation_factor)
        t=correct_temperature(v=t,vo=vo)
    else:
        t=inputT.copy()
    
    if inputCompression==None:
        compression=decimate(data[start_point:end_point,pot1]+data[start_point:end_point,pot2].astype("d"),decimation_factor)
    else:
        compression=inputCompression.copy()
    
    summary_stats=[]
    maxr=0
    for smooth in np.arange(smooth_min,smooth_max,smooth_step):
        t_smooth=smooth_data(t,smooth)
        poly=np.polyfit(t_smooth[smooth:],compression[smooth:],3)
        x=t_smooth[smooth:]
        r=np.corrcoef(poly[0]*x**3+poly[1]*x**2+poly[2]*x+poly[3],compression[smooth:])[0,1]
        if r>maxr:
            maxr=r
            summary_stats=[smooth,r,poly]

    smooth=summary_stats[0]
    final_t=smooth_data(t,smooth)
    fitt=np.zeros(final_t.shape)
    poly=summary_stats[2]
    for i in range(len(poly)):
        fitt+=poly[len(poly)-i-1]*final_t**i
    return summary_stats,compression[smooth:],t[smooth:],final_t[smooth:],fitt[smooth:]

def main():
        
    """docstring for main"""
    print("Loading data")
    data=swim_io.read_nc("full_interception.nc").data
    if master_end_point==None:
        end_point=data.shape[0]
    elif type(master_end_point)==float:
        if master_end_point<1:
            end_point=int(data.shape[0]*master_end_point)
        else:
            end_point=master_end_point
    else:
        end_point=master_end_point
    print("Fitting T")
    stats,compression,rawt,tsmooth,fitt=optimal_t_fit(data)
    print("Optimal smoothing:"+str(stats[0]))
    r=str(stats[1])[:4]
    print("r :"+r)
    
    print("Applying to data")
    v=decimate(data[:,8].astype("d"),decimation_factor)
    vo=decimate(data[:,-1].astype("d"),decimation_factor)
    rawt=correct_temperature(vv=v/vo)
    rawcompression=decimate(data[:,pot1]+data[:,pot2].astype("d"),decimation_factor)
    full_smootht=smooth_data(rawt,stats[0])
    full_fitt=np.zeros(full_smootht.shape)
    poly=stats[2]
    for i in range(len(poly)):
        full_fitt+=poly[len(poly)-i-1]*full_smootht**i
    
    print("Calculating dates")
    startdate=datetime.datetime(2013,10,3,10,0)
    enddate=datetime.datetime(2014,3,5,13,0,0)
    startdate=datetime.datetime(2014,4,28,11,20,0)
    enddate=datetime.datetime(2014,5,16,12,27,0)
    dt=(enddate-startdate)/len(rawt)
    dates=[startdate+dt*i for i in range(len(rawt))]
    subdates=dates[start_point/decimation_factor:end_point/decimation_factor][stats[0]:-1]
    
    print("Plotting...")
    plt.figure(figsize=(14,6),dpi=100);
    plt.plot(dates,rawcompression,label="Compression")
    plt.plot(dates,full_fitt,label="Fitted Temp.")
    plt.plot(dates,rawcompression-full_fitt+rawcompression[len(rawcompression)/2],label="T.Corr. Comp.")
    plt.plot([dates[0],dates[1]],[-100,-100],label="Temperature",color="grey")
    plt.ylabel("Compression")
    if tree==1:
        plt.ylim(25000,27000)
    if tree==2:
        plt.ylim(27000,35000)
    if tree==3:
        plt.ylim(27000,30000)
    # plt.ylim(20000,30000)
    plt.legend()
    plt.title("Tree:{}  Potentiometers:{} & {}  Temp.fit:r={}".format(tree,pot1,pot2,r))
    # plotting temperature on a second axis
    temp_axis=plt.twinx()
    temp_axis.plot(dates,rawt,color="grey")
    temp_axis.set_ylim(10,-50)
    temp_axis.set_ylabel("Temperature [C]")
    date_index=start_point/decimation_factor
    temp_axis.plot([dates[date_index],dates[date_index]],[10,-50],"--",color="black")
    date_index=end_point/decimation_factor
    if date_index>=len(dates):date_index=-1
    temp_axis.plot([dates[date_index],dates[date_index]],[10,-50],"--",color="black")
    temp_axis.plot([dates[0],dates[-1]],[0,0],color="red")
    
    plt.xlim(dates[int(len(dates)*0.05)],dates[int(len(dates)*0.95)])
    plt.savefig("temp-v-compression_tree{}.png".format(tree))

if __name__ == '__main__':
    
    main()