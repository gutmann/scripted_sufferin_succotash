#!/usr/bin/env python
# encoding: utf-8
"""
process_snow_depth.py

Created by Ethan Gutmann on 2011-03-16.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""

#Standard
import sys
import getopt
import pdb
#Library
import numpy as np
from numpy import sqrt,where,array
import numpy.ma as ma
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as font_prop
#3rd party
import julday
#user
import load_data
import date_fun



help_message = '''
Processes snow depth data, calibrating and removing temperature effects

USAGE:
    process_snow_depth datafile keyfile calfile [badfile] [--inches]
    
    datafile = input data file name
    keyfile = main metadata file
            Column formatted text (with arbitrary non-numeric header)
                Column 1: Instrument number (must be a number, not Judd 7-2 or even 7-2)
                Column 2: Judd distance [cm] (if distance not supplied enter -99)
                Column 3: Judd Depth [cm] (if depth not supplied, enter -99)
                Column 4: Judd time [ms] (if time not supplied, enter -99)
                Column 5: Judd Temperature [deg. C] (if supplied, data will be corrected)
                Column 6: Sensor location, 1=open 2=sub-canopy (but could be used for any two subsets)
    calfile = main calibration data file
            Column formatted text (with arbitrary non-numeric header)
                First line = identifies each column, 
                    negative numbers for dates:
                        year=-1,month=-2,...minute=-5
                    positive numbers for instrument numbers (corresponding to the key file instrument numbers)
                    Example: -1 -2 -3 -4 -5 1 2 3 4 5 6
                subsequent lines = data
                    Example: 2011  1  23  18  0      99.3    75.0    86.7    82.3    101.7   94.3
    badfile = optional description of periods of bad data
            Column formatted text (with arbitrary non-numeric header)
            columns are : Instrument number,start Year,month,day,hour,minute end year,month,day,hour,minute
            Example:
                Stn     Start (Y,M,D,h,m)   End-Bad (Y,M,D,h,m)
                4       2011 2 21 2 0       2011 5 8 19 0
    -i,--inches = set this flag if the input data are in inches 
            data will be converted to cm, calibration data should always be in cm
    -f,--Farenheight = set this flag if the input air temperatures are in Farenheight 
    -c,--chimney = set this flag if the input data are from Chimney park in 2010-2011, 
            a fix for the flip-flop of time and distance is corrected

NOTE: the order of chimney and farenheight options seems to be important!           
EXAMPLE:
    process_snow_depth.py --chimney --Farenheight grey_proc.txt grey_keys.txt grey_cal.txt
'''
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def apply_calib(stimes,snow,ctimes,cal):
    
    npoints=len(cal)
    offsets=np.zeros(npoints)
    for i in range(npoints):
        curpos=where(abs(stimes - ctimes[i]) < 1E-3)[0]
        offsets[i]=cal[i]-snow[curpos]

        if i == 0 :
            snow[0:curpos]+=offsets[i]
        else:
            n=curpos-lastpos
            delta=(offsets[i]-offsets[i-1])/float(n)
            if (delta==0):
                curoffset=offsets[i]
            else:
                curoffset=(np.arange(offsets[i-1],offsets[i],delta))[0:n]
            snow[lastpos:curpos]+=curoffset

        lastpos=curpos
        
    if curpos < len(snow):
        snow[curpos:]+=offsets[-1]
    
    return snow
        

def simple_cleanup(data,mask):
    """
    Simply removes all values above or below two thresholds and substitutes with the previous good value
    
    """
    max_condition=140.0
    min_condition=-10.0
    i=0
    if (min_condition > data[i] > max_condition):
        tmp=where((data >= min_condition) | (data<=max_condition))[0]
        data[0]=data[tmp[0]]
        mask[0]=0
    
    
    for i in np.arange(len(data) - 1) + 1:
        if (min_condition > data[i] > max_condition):
            data[i]=data[i-1]
            mask[i]=0
        if mask[i-1]==1:
             if (np.abs(data[i]-data[i-1]) > 15):
                data[i]=data[i-1]
                mask[i]=0
    
    morebad=where(~np.isfinite(data))
    data[morebad] = -99
    mask[morebad] = 0
    return data



def fix_badsnow(snow,mask):
    """
    find regions marked "bad" and fill them in using a linear regression
    """
    # from numpy.linalg import lstsq
    from ols import ols
    
    # find times when all of the snow data are "good"
    # allgood=where(np.min(snow,axis=1) >= 0)[0]
    
    # loop over all snow columns fixing bad data
    for i in range(len(snow[0,:])):
        # find bad data
        bad=where(mask[:,i] == 0)[0]
        # if there is any bad data fix it
        if len(bad) >0:
            # cursnow is data to be fixed
            cursnow=snow[:,i]
            # othersnow is data to use to fix it
            othersnow=np.delete(snow,i,1)
            othermask=np.delete(mask,i,1)
            
            for j in range(len(bad)):
                otherbadsnows=where(othermask[bad[j],:]==0)[0]
                if len(otherbadsnows)<len(othersnow[0,:]):
                    if len(otherbadsnows)>0:
                        useothersnow=np.delete(othersnow,otherbadsnows,1)
                        useothermask=np.delete(othermask,otherbadsnows,1)
                    else: 
                        useothersnow=othersnow
                        useothermask=othermask
                    
                    # print("available data="+str(useothersnow.shape))
                    # find times when all of the other snow data are "good"
                    allgood=where((np.min(useothermask,axis=1) > 0) & (cursnow>0))[0]
                    # print("good data points="+str(len(allgood)))
                    # perform linear regression on the remaining data using the last and next 10 days of data
                    usepoints= (24*4*10l)
                    if len(allgood)>24:
                        nearestpoints=where(np.abs(allgood-bad[j])<usepoints)[0]
                        if len(nearestpoints)<400:
                            nearestpoints=where(np.abs(allgood-bad[j])<usepoints*5)[0]
                            # print("using more points")
                        if len(nearestpoints)>=500:
                            line=ols(cursnow[allgood[nearestpoints]],useothersnow[allgood[nearestpoints],:])
                            snow[bad[j],i]=line.b[0]+np.sum(useothersnow[bad[j],:]*line.b[1:])
                            # print("Enough points: Station "+str(i)+"  time: "+str(bad[j]))
                        else:
                            # print("Station "+str(i)+"  time: "+str(bad[j])+"  oldval:"+str(snow[bad[j],i])+"  Newval:"+str(snow[bad[j]-1,i]))
                            snow[bad[j],i]=snow[bad[j]-1,i]
                    else:
                        snow[bad[j],i]=snow[bad[j]-1,i]
            
            
            
            
            # figure out which columns have bad data that overlap with cursnow's bad data
            # otherbadsnows=where(np.min(othersnow[bad,:], axis=0) < 0)[0]
            # if there are columns with bad data in this period, remove them
            # if len(otherbadsnows)>0:
            #   othersnow=np.delete(othersnow,otherbadsnows,1)
                
            # perform linear regression on the remaining data using only the last 10 days of data
            # usepoints= (24*4*10l)
            # usepoints=len(allgood)-usepoints
            # usepoints=0
            
            # import pdb
            # pdb.set_trace()
            
            # line=ols(cursnow[allgood[usepoints:]],othersnow[allgood[usepoints:],:])
            # start with the offset
            # print(line.b)
            # snow[bad,i]=line.b[0]
            # loop over 
            # for thatsnow in range(len(othersnow[0,:])):
            #   snow[bad,i]+=othersnow[bad,thatsnow]*line.b[thatsnow+1]

        
# pass a moving median filter over the entire timeseries using a window size specified below
#  initial windowsize = 16 (w/15min dt = 4hr)
def smooth_snow(snow):
    sz=snow.shape
    
    windowsize=4
    halfwin=windowsize/2
    for i in np.arange(sz[0]-windowsize-1)+halfwin:
        for j in range(sz[1]):
            snow[i,j]=np.median(snow[i-halfwin:i+halfwin+1,j])


def make_times(dates):
    output=np.zeros((len(dates[:,0])))
    for i in range(len(dates[:,0])):
        output[i]=julday.mjul_day(dates[i,0],dates[i,1],dates[i,2],dates[i,3],dates[i,4],0)
        
    return output

def remove_baddays(snow,dtimes,badtimes,mask):
    
    for i in range(len(badtimes[:,0])):
        temp=where((badtimes[i,0]<dtimes)&(dtimes<badtimes[i,1]))[0]
        if len(temp) >0:
            snow[temp]=-99
            mask[temp]=0


def main(datafile, keyfile, calfile,badfile,inches=False,Farenheight=False,Chimney_2011=False):
    data=load_data.cols(datafile,dtype='d')
    keys=load_data.cols(keyfile)
    if badfile!=None:
        baddata=load_data.cols(badfile)
        badtimes=np.array((make_times(baddata[:,1:6]),make_times(baddata[:,6:]))).T
    
    (vtimes,val)=load_data.cols_date(calfile,year=0,month=1,day=2,hour=3,minute=4)
    
    vtimes=vtimes[1:]
    # dtimes=date_fun.excel2mjd(data[:,0])
    dtimes=data[:,0]
    dates=date_fun.mjd2datetime(dtimes,roundseconds=True)
    
    allsnow=np.zeros((len(dtimes),len(keys[:,0])))
    allmask=np.zeros((len(dtimes),len(keys[:,0])))
    fig=plt.figure()
    
    for i in range(len(keys[:,0])):
        if keys[i,1]>0:
            snow=data[:,keys[i,1]]
        if keys[i,3]>0: 
            snow=data[:,keys[i,3]]/6.02
        if keys[i,2]>0:
            snow=100-data[:,keys[i,2]]
        if keys[i,4]>0:
            if Farenheight:
                print("Working in Farenheight")
                airt=(data[:,keys[i,4]]-32)/1.8+273.15
            else: 
                print("we shouldn't get here for chimney park")
                airt=data[:,keys[i,4]]+273.15
            tmp=np.where(airt<0)
            if len(tmp[0])>0:
                airt[tmp]=9999
                snow[tmp]=-9999
            if not Chimney_2011:
                print("we shouldn't get here for chimney park")
                snow*=sqrt(airt/273.15)
        if inches and not Chimney_2011:
            print("we shouldn't get here for chimney park")
            print("Working in Inches")
            snow*=2.54

        if Chimney_2011:
            print("Correcting Chimney Park")
            snow*=25.4 # "snow" is actually time that was mistakenly converted to inches from mm on the datalogger.
            snow*=0.3314/2 #0.3314 cm/microsecond = speed of sound in air at 0C (air temperature correction follows)
            snow*=sqrt(airt/273.15)
    
        # dists=snow.copy()
        mask=np.ones(len(snow))
        tmp=where(snow<50)
        snow=200-snow
        if len(tmp[0])>0:
            snow[tmp]=-99
            mask[tmp]=0
        
        thesekeys=where(val[0,:] == keys[i,0])
        curval=val[1:,thesekeys]
        curval.shape=(len(curval))
        good_val=where(curval >= 0)
        # this was used for 2010 data from one site where it needed to be cleaned up before calibration
        # if i==1:snow=simple_cleanup(snow,mask) 
        
        snow=apply_calib(dtimes, snow,vtimes[good_val[0]],curval[good_val[0]])
        snow=simple_cleanup(snow,mask)
        tmp=where(snow<-10)
        if len(tmp[0])>0:
            snow[tmp]=np.nan
            plt.plot(dates,snow)
            snow[tmp]=-99
            mask[tmp]=0
        else:
            plt.plot(dates,snow)
        plt.xlabel("Date")
        plt.ylabel("Snow Depth (cm)")
        fig.autofmt_xdate()
        fig.savefig('snow_depth_'+str(int(keys[i,0]))+'.pdf')
        fig.clf()
        
        if badfile!=None:
            baddays=where(baddata[:,0] == keys[i,0])
            if len(baddays) >0 :
                for thisday in baddays:
                    remove_baddays(snow,dtimes,badtimes[thisday,:],mask)
        
        allmask[:,i]=mask
        allsnow[:,i]=snow
    
    # allsnow[:,-1]=0
    # allmask[:,-1]=0
    tmp=where(allmask==0)
    if len(tmp[0])>0:allsnow[tmp]=-99
    # fix_badsnow(allsnow[:,:-1],allmask[:,:-1])
    fix_badsnow(allsnow,allmask)
    allsnow[where(allsnow<0)]=0

    # tmp=where(dtimes >julday.mjul_day(2011,6,21,0,0,0))[0]
    # allsnow[tmp[0]:,:]=0
    
    morebaddata=where(~np.isfinite(allsnow))
    if len(morebaddata[0])>0:
        allsnow[morebaddata]=-99
    smooth_snow(allsnow)

    for i in range(len(allsnow[0,:])):
        plt.plot(dates,allsnow[:,i])
        plt.xlabel("Date")
        plt.ylabel("Snow Depth (cm)")
        fig.autofmt_xdate()
        fig.savefig('fixed_snow_depth_'+str(int(keys[i,0]))+'.pdf')
        fig.clf()
    
    
    for i in range(len(keys[:,0])):
        with open('snow_depth_'+str(int(keys[i,0]))+'.txt', 'w') as f:
            f.write('   Date      Time    Snow Depth    QCflag\n')
            f.write('yyyy-mm-dd hh:mm:ss     (cm)    (0=bad,1=good)\n')
            [f.write(str(dates[j])+'      '+str(round(allsnow[j,i]*100)/100.0)+'          '+
            str(int(allmask[j,i]))+'\n') for j in range(len(snow))]
    
    with open('snow_depth_all.csv','w') as f:
        f.write('   Date      Time  ')
        for i in range(len(keys[:,0])):
            f.write(', Snow Depth,   QCflag     ')
        f.write('\n')
        
        f.write('yyyy-mm-dd hh:mm:ss')
        for i in range(len(keys[:,0])):
            f.write(',    (cm),   (0=bad 1=good)')
        f.write('\n')

        for j in range(len(snow)):
            f.write(str(dates[j])+',')
            for i in range(len(keys[:,0])):
                f.write('     '+str(round(100*allsnow[j,i])/100.0)+',         '+
                        str(int(allmask[j,i]))+',        ')
            f.write('\n')
    

    
    color_vals=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(0.5,1,0),(0.5,0,1),(0,0.5,1),(1,0.5,0),(1,0,0.5),(0,1,0.5)]
    for i in range(len(keys[:,0])):
        plt.plot(dates,allsnow[:,i],label='Judd-'+str(int(keys[i,0])),lw=0.5,color=color_vals[i])
        
    plt.plot(dates, np.mean(allsnow,axis=1),lw=3.0,label='Mean',color='k')
    # import pdb; pdb.set_trace()
    i=len(keys[:,0])+1
    plt.plot(dates, np.mean(allsnow[:,where(keys[:,5] ==1)[0]],axis=1),lw=2.0,ls='dashed',label='Open',color=color_vals[i])
    i=i+1
    plt.plot(dates, np.mean(allsnow[:,where(keys[:,5] ==2)[0]],axis=1),lw=2.0,ls='dashed',label='Sub-Canopy',color=color_vals[i])
    plt.legend(ncol=4,prop=font_prop(size=11),loc=2)
    plt.xlabel("Date")
    plt.ylabel("Snow Depth (cm)")
    fig.autofmt_xdate()
    plt.ylim(0,200)
    plt.show()
    fig.savefig('plot_summary.pdf')

def go(f1,f2,f3,f4):
    main(f1,f2,f3,f4)

def cli_main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hoicf:v", ["help", "output=","inches","chimney","Farenheight"])
        except getopt.error, msg:
            raise Usage(msg)
        
        Farenheight=False
        Chimney_2011=False
        inches=False
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
            if option in ("-i", "--inches"):
                inches=True
            if option in ("-f", "--Farenheight"):
                Farenheight=True
            if option in ("-c", "--chimney"):
                Chimney_2011=True
            print(option)
        print(opts)
                
        if len(args)<3:
            raise Usage(help_message)

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2

    print(Farenheight)
    if len(args)==3:
        main(args[0],args[1],args[2],None,inches=inches,Farenheight=Farenheight,Chimney_2011=Chimney_2011)
    else:
        main(args[0],args[1],args[2],args[3],inches=inches,Farenheight=Farenheight,Chimney_2011=Chimney_2011)


if __name__ == "__main__":
    sys.exit(cli_main())
