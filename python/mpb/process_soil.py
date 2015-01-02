#!/usr/bin/env python
# encoding: utf-8
"""
process_soil.py

Created by Ethan Gutmann on 2011-04-06.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""
#Library
import numpy as np
from numpy import sqrt,where,array
import numpy.ma as ma
import matplotlib
matplotlib.use('PDF') #should only be performed if not in an interactive session
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as font_prop
#3rd party
import julday #shouldn't use this anymore, just use the date_fun functions
#user
import load_data
import date_fun

def simple_cleanup(data):
    """
    Simply removes all values above or below two thresholds
    """
    max_condition=0.55
    min_condition=0.0
    i=0
    # if (min_condition > data[i] > max_condition):
    #   tmp=where((data >= min_condition) | (data<=max_condition))[0]
    #   data[0]=data[tmp[0]]
    # 
    # for i in np.arange(len(data) - 1) + 1:
    #   if (min_condition > data[i] > max_condition):
    #       data[i]=data[i-1] 
        
    data[where(data > max_condition)] = -9999
    data[where(data < min_condition)] = -9999
    
    return data

def cleanup(smc):
    '''
    clean up soil moisture spikes
    
    search through for points in which soil moisture changes by >0.01 
    immediately before *and* after the current value
    replace those spikes with the mean of the surrounding values
    '''
    # for i in np.arange(smc.shape[1]):
    for j in np.arange(smc.shape[0]-2)+1:
        if ((abs(smc[j]-smc[j-1])>0.1) & 
            (abs(smc[j]-smc[j+1])>0.1)):
            smc[j]=(smc[j-1]+smc[j+1])/2
            
    max_condition=0.5
    min_condition=0.0
    smc[where(smc > max_condition)] = -9999
    smc[where(smc < min_condition)] = -9999
    


    
    
def simple_fillsmc(smc):
    '''
    fill missing soil moisture values
    
    search through the soil moisture array for missing (<0) values
    when one is found, find the next valid value (>0) and interpolate linearly
    between valid values. 
    '''
    for i in np.arange(smc.shape[1]):
        j=0
        while j<smc.shape[0]-1:
            current=smc[j,i]
            if current<0:
                startpoint=j-1
                while (current<0) and (j<smc.shape[0]-1):
                    j+=1
                    current=smc[j,i]
                endpoint=j+1
                if (current>0) and (startpoint>=0):
                    startsmc=smc[startpoint,i]
                    endsmc=smc[endpoint,i]
                    
                    smc[startpoint+1:endpoint,i]=startsmc+(endsmc-startsmc)*np.arange(endpoint-startpoint-1)
                
            else:
                pass
            j+=1
    


def fix_badsmc(smc):
    # from numpy.linalg import lstsq
    from ols import ols
    smc_standard=smc.copy()
    # find times when all of the smc data are "good"
    allgood=where(np.min(smc,axis=1) >= 0)[0]
    while len(allgood)<(24*4*15):
        sumall=np.sum(smc_standard,axis=0)
        worse=np.min(sumall)
        betterprobes=np.where(sumall>sumall[worse])
        if len(betterprobes[0])<2:
            return
        smc_standard=smc_standard[:,betterprobes]
        allgood=where(np.min(smc,axis=1) >= 0)[0]
        
    
    # loop over all smc columns fixing bad data
    for i in range(len(smc[0,:])):
        # find bad data
        bad=where(smc[:,i] < 0)[0]
        # if there is any bad data fix it
        if len(bad) >0:
            print(len(allgood))
            # cursmc is data to be fixed
            cursmc=smc[:,i]
            # othersmc is data to use to fix it
            othersmc=np.delete(smc,i,1)
            # figure out which columns have bad data that overlap with cursmc's bad data
            otherbadsmcs=where(np.min(othersmc[bad,:], axis=0) < 0)[0]
            # if there are columns with bad data in this period, remove them
            if len(otherbadsmcs)>0:
                othersmc=np.delete(othersmc,otherbadsmcs,1)
                
            # perform linear regression on the remaining data using only the last 15days? of data
            # usepoints= (-1*24*4*15)
            usepoints=0
            if len(allgood)>50:
                print(len(allgood))
                line=ols(cursmc[allgood[usepoints:]],othersmc[allgood[usepoints:],:])
                # start with the offset
                smc[bad,i]=line.b[0]
                # loop over 
                for thatsmc in range(len(othersmc[0,:])):
                    smc[bad,i]+=othersmc[bad,thatsmc]*line.b[thatsmc+1]

        
def smooth_smc(data):
    """
    Pass a moving median filter over the entire timeseries using a window size specified below
    initial windowsize = 16 (w/15min dt = 4hr)
    windowsize=4 => 1hr
    add one to the end step so it will be centered on the value being filtered
    """
    sz=data.shape
    
    windowsize=4
    halfwin=windowsize/2
    for i in np.arange(sz[0]-windowsize-1)+halfwin:
        for j in range(sz[1]):
            data[i,j]=np.median(data[i-halfwin:i+halfwin+1,j])


def make_times(dates):
    '''convert an array (n x [year,month,day,hour,minute]) to modified Julian days'''
    output=np.zeros((len(dates[:,0])))
    for i in range(len(dates[:,0])):
        output[i]=julday.mjul_day(dates[i,0],dates[i,1],dates[i,2],dates[i,3],dates[i,4],0)
    return output


def remove_baddays(smc,dtimes,badtimes):
    '''find dates that fall within a "badtime" and remove soil moisture data for that period'''
    for i in range(len(badtimes[:,0])):
        temp=where((badtimes[i,0]<dtimes)&(dtimes<badtimes[i,1]))[0]
        if len(temp)>0:
            # print("removing:"+str(date_fun.mjd2datetime(badtimes[i,0]))+" -- "+str(date_fun.mjd2datetime(badtimes[i,1])))
            smc[temp]=-9999

def remove_frozen(smc,tsoil):
    '''find times when the soil is frozen and remove the soil moisture data from that peroiod'''
    temp=where(tsoil<0.1)[0]
    if len(temp) >0:smc[temp]=-9999

def make_plots(allsmc,dates,keys):
    ''' make output pdf plots'''
    fig=plt.figure()
    color_vals=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),
                (0.5,1,0),(0.5,0,1),(0,0.5,1),(1,0.5,0),(1,0,0.5),(0,1,0.5),
                (0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5)]
    
    # plot all the raw data
    for i in range(len(keys[:,0])):
        plt.plot(dates,allsmc[:,i],label='SMC-'+str(int(keys[i,0])),lw=0.5,color=color_vals[i])
    # now plot the means for each soil level
    plt.plot(dates, np.mean(allsmc,axis=1),lw=3.0,label='Mean',color='k')
    i=len(keys[:,0])+1
    plt.plot(dates, np.mean(allsmc[:,where(keys[:,4] ==1)[0]],axis=1),lw=2.0,
            ls='dashed',label='10cm',color=color_vals[i])
    i=i+1
    plt.plot(dates, np.mean(allsmc[:,where(keys[:,4] ==2)[0]],axis=1),lw=2.0,
            ls='dashed',label='30cm',color=color_vals[i])
    i=i+1
    plt.plot(dates, np.mean(allsmc[:,where(keys[:,4] ==3)[0]],axis=1),lw=2.0,
            ls='dashed',label='60cm',color=color_vals[i])
    # finally, make a legend, add axis-labels, and save the figure
    plt.legend(ncol=5,prop=font_prop(size=11),loc=2)
    plt.xlabel("Date")
    plt.ylabel("Soil Moisture [$cm^3/cm^3$]")
    fig.autofmt_xdate()
    plt.ylim(0.00,0.55)
    fig.savefig('soil_moisture_all.pdf')
    
    # plot just the means so they can be seen without the clutter of the raw data
    fig.clf()
    i=len(keys[:,0])+1
    plt.plot(dates, np.mean(allsmc[:,where(keys[:,4] ==1)[0]],axis=1),lw=2.0,
            ls='dashed',label='10cm',color=color_vals[i])
    i=i+1
    plt.plot(dates, np.mean(allsmc[:,where(keys[:,4] ==2)[0]],axis=1),lw=2.0,
            ls='dashed',label='30cm',color=color_vals[i])
    i=i+1
    plt.plot(dates, np.mean(allsmc[:,where(keys[:,4] ==3)[0]],axis=1),lw=2.0,
            ls='dashed',label='60cm',color=color_vals[i])
    # again, make the legend, axis-labels, and save the figure
    plt.legend(ncol=5,prop=font_prop(size=11),loc=2)
    plt.xlabel("Date")
    plt.ylabel("Soil Moisture [$cm^3/cm^3$]")
    fig.autofmt_xdate()
    plt.ylim(0.00,0.55)
    fig.savefig('soil_moisture_means.pdf')
    # return interactive plot control(?)
    plt.show()


def main(datafile, keyfile, badfile=None,topp=False):
    '''
    cleanup soil moisture data, fill bad data, add QC column
    
    to edit: 
        Temperature = [C] or [K]? (currently [C] : line 265 )
        date = excel or mjd or ... (currently excel : line 258)
    
    process_soil takes a datafile and a keyfile as input (and optionally a badfile)
    The data file should be a column formatted file with 
        Column 1 = Modified Julian Day or excel date (must modify code)
        Column 2-n = soil moisture and temperature (others ignored)
    The key file should be a column formatted file too with:
        Column 1: Instrument number (must be a number, not Judd 7-2 or even 7-2)
        Column 2: Soil Moisture Probe Column number (0 based)
        Column 3: Soil Temperature Probe Column number (0 based)
        Column 4: Soil Moisture - Temperature calibration (cm3/cm3/C)
        Column 5: Sensor group, 1=10cm 2=30cm 3=60cm (but could be used for any subsets)
        
    badfile is also column formatted specifying beginning and ending dates of bad data periods
    
    Soil moisture data is first roughly QAQCed by removing spikes, and data values outside
    predifined thresholds.  
    '''
    
    # load the data files or fail
    try:
        data=load_data.cols(datafile,dtype='d')
        keys=load_data.cols(keyfile)
    except IOError: 
        print("Badly formed input data file")
        return
    # data[:,3]=-9999
    # remove any rows in which the date is negative
    tmp=where(data[:,0]>0)[0]
    if len(tmp)<10: 
        print("No valid dates in file")
        return
    data=data[tmp,:]
    
    # create datetime objects for use plotting and printing "pretty"
    dtimes=data[:,0] #if dates are in mjd this is all you need
    # dtimes=date_fun.excel2mjd(dtimes) #if dates are in excel format, convert them to mjd
    # dtimes=make_times(data[:,0:5]) #if dates are in year,month,day,hour,minute convert them to mjd
    dates=date_fun.mjd2datetime(dtimes,roundseconds=True)
    # topp=True
    
    # set up an array to hold all of the soil moisture data
    allsmc=np.zeros((len(dtimes),len(keys[:,0])))
    # loop over all soil moisture columns defined in the keyfile
    for i in range(len(keys[:,0])):
        # grab the current soil moisture and temperature data 
        smc=data[:,keys[i,1]]
        if topp==True:
            print("Applying Topp calibration")
            smc=-5.3e-2+2.92e-2*smc-5.5e-4*smc**2+4.3e-6*smc**3
        # data[:,keys[i,2]]-=273.15
        tsoil=data[:,keys[i,2]]
        if np.median(tsoil)>200:
            tsoil-=273.15
        # find where the soil temperature is good and apply a temperature correction
        # to the soil moisture data
        tmp=where(tsoil >-100)[0]
        smc[tmp]=smc[tmp]-((tsoil[tmp]-10)*keys[i,3])
        # perform some minimal cleanup on the soil data (remove spikes and out of bounds)
        cleanup(smc)
        
        remove_frozen(smc,tsoil)
        
        if badfile != None:
            baddata=load_data.cols(badfile)
            if baddata.shape[1]>3:
                badtimes=np.array([make_times(baddata[:,1:6]),make_times(baddata[:,6:])]).T
            else:
                badtimes=baddata[:,1:]
            baddays=where(baddata[:,0] == keys[i,0])
            if len(baddays) >0 :
                for thisday in baddays:
                    remove_baddays(smc,dtimes,badtimes[thisday,:])
        
        allsmc[:,i]=smc
        
    # create a QC mask before fixing the bad soil moisture values
    # allsmc[where(allsmc<0.0005)]=-9999
    allmask=np.ones(allsmc.shape)
    allmask[where(allsmc<0)]=0
    # now fix bad soil moisture values
    fix_badsmc(allsmc)
    # use a median filter to remove small spikes. 
    # smooth_smc(allsmc)
    
    # make all PDF plots
    make_plots(allsmc,dates,keys)
    # write output file (could be moved to a subroutine)
    with open(datafile+'_all.csv','w') as f:
        f.write('   Date      Time  '+',     Soil Temp.,  Moisture, QCflag     '*len(keys[:,0])+'\n')
        f.write('yyyy-mm-dd hh:mm:ss'+',         (C),   (   cm), (0=bad 1=good)'*len(keys[:,0])+'\n')
        
        for j in range(len(allsmc[:,0])):
            f.write(str(dates[j])+',')
            for i in range(len(keys[:,0])):
                f.write('     %8.3f,  %7.3f,   %3i,         ' %(data[j,keys[i,2]],allsmc[j,i],allmask[j,i]))
            f.write('\n')

# print commandline documentation
def print_doc():
    print("""
process_soil.py datafile keyfile [badfile]
    python script to process a set of column formatted soil moisture data
    datafile should be space or tab separated columns of data
    keyfile should be space or tab separated columns of metadata in the form
        identifier SMC_column Temp_column smc-Ts_calibration
    optionally: badfile is a space or tabe separated 
        
for keyfile: 
    identifier is an arbitrary number identifying each soil moisture data set uniquely
    smc_column is the column number (zero based) for this soil moisture data set
    Temp_column is the column number (zero based) for the corresponding soil temperature
    smc-Ts_calibration is the soil moisture -temperature calibration relationship (cm^3/cm^3)/DegC (typically ~0.001-0.0003)
    """)

# if called from a UNIX Shell, parse args and call main
if __name__ == '__main__':
    import sys
    
    if len(sys.argv) >3:
        main(sys.argv[1],sys.argv[2],sys.argv[3])
    elif len(sys.argv) >2:
        main(sys.argv[1],sys.argv[2])
    else: print_doc()

