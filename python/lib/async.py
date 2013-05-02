from __future__ import print_function
import sys
import numpy as np
from scipy import stats
from bunch import Bunch
# from nc_reader import NC_Reader
# import mygis
narrfilename='/Volumes/Data2/usbr/NARR/flxnc/HGT_2006010100.nc'
'''Default filename to use for low res geographic information'''
wrffilename="/Volumes/G-SAFE/usbr/wrf4km-sfc-output/wrfout_d01_2008-01-01_00:00:00.nc"
'''Default filename to use for high res geographic information'''
narrdatafiles="/Volumes/Data2/usbr/NARR/precip/merged_AWIP32.2008*"
'''Default search path to use for low res data sets'''
wrfdatafiles=["/Volumes/G-SAFE/usbr/wrf4km_daily_precip/NARR_04km_OCT2007-SEP2008.nc",
              "/Volumes/G-SAFE/usbr/wrf4km_daily_precip/NARR_04km_OCT2008-DEC2008.nc"]
'''Default files to use for high res data'''

def read_wrf():
    '''Read in WRF precip data'''
    pass
    # # read data from two files, the first is a full Water Year, 
    # # the second completes the calendar year of the first
    # d1=mygis.read_nc(wrfdatafiles[0],var="RAINNC").data
    # d2=mygis.read_nc(wrfdatafiles[1],var="RAINNC").data
    # # nday shift = the number of days to strip from the beginning of the first file
    # nday_shift=d2.shape[0] # (31+30+31)
    # # compute the first day from accumulated values
    # day1=d1[nday_shift,:,:]-d1[nday_shift-1,:,:]
    # # strip off / shift data by n_days
    # d1[:-nday_shift,:,:]=d1[nday_shift:,:,:]
    # # tack on the second files data to the end
    # d1[-nday_shift:,:,:]=d2
    # # compute the difference from accumulated precip for the rest of the data
    # d1[1:,:,:]=np.diff(d1,axis=0)
    # # finally restore the first day computed previously
    # d1[0,:,:]=day1
    # return d1

def read_narr():
    '''Read in NARR precip dataset from N files into an array(N,y,x) matching a hires file'''
    # make an NC_Reader object that can loop over all files and geographically match them
    # to the wrf input data file
    pass
    # narr=NC_Reader(narrdatafiles,geoin_file=narrfilename,geomatch_file=wrffilename,
    #                 readvars=list(['PRATE_221_SFC']),nn=False, bilin=True)
    # # read the first (and only) variable ([0]) from the first data file
    # narrd1=narr.next()[0]
    # # create an output array to handle the data
    # narroutputdata=np.empty((len(narr._filenames),narrd1.shape[0],narrd1.shape[1]),dtype=np.float32)
    # # insert the first day
    # narroutputdata[0,:,:]=narrd1*3*60*60
    # # could probably do this as an enumerate(narr) ?
    # # i is the index into the output dataset
    # i=1
    # for d in narr:
    #     # loop through the narr dataset storing each value in the outputdata. 
    #     narroutputdata[i,:,:]=d[0]*3*60*60 #convert mm/s to mm over a three hour period
    #     i+=1
    # # now strip the data down from 3hrly to daily data (3hr = 8x per day)
    # daily=np.empty((len(narr._filenames)/8,narrd1.shape[0],narrd1.shape[1]),dtype=np.float32)
    # for i in range(daily.shape[0]):
    #     daily[i]=narroutputdata[i*8:(i+1)*8,:,:].sum(axis=0)
    #     
    # return daily

def make_xy(inx,iny,precip=True,n=300.0):
    '''Make x,y into <= n [300] roughly evenly spaced values'''
    # sort data so we can perform an asynchronous regression
    n=np.float(n)
    inx.sort()
    iny.sort()
    # convert to log for precip values
    if precip:
        # set minimum values so the log doesn't explode
        inx[inx<1E-3]=1E-3
        iny[iny<1E-3]=1E-3
        inx=np.log(inx)
        iny=np.log(iny)
    # find min and max values to define range and step size for evenly spaced output
    maxx=np.max([inx[-1],iny[-1]])
    minx=np.min([inx[0],iny[0]])
    fullrange=maxx-minx
    step=(fullrange)/(n)
    # create output variables
    outy=np.empty(n,dtype=np.float32)
    outx=np.empty(n,dtype=np.float32)# np.arange(minx,maxx,step)+step/2
    for i in range(int(n)):
        # set up current upper and lower bounds
        lower=(i/n)*fullrange+minx
        upper=((i+1)/n)*fullrange+minx
        # find mean of all values with X values in that range
        curvals=np.where((inx<=upper)&(inx>lower))
        if len(curvals[0])>0:
            outy[i]=np.mean(iny[curvals])
            outx[i]=np.mean(inx[np.where((inx<=upper)&(inx>lower))])
        else:
            outy[i]=np.nan
            outx[i]=np.nan
    # remove bad values (where no valid data existed)
    tmp=np.where(np.isfinite(outy))
    return(outx[tmp],outy[tmp])
        
def test_plot(x,y,linfits,i,j):
    import time
    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(x,y,".")
    for i in range(linfits.shape[1]):
        x1,x2=linfits[0,i],linfits[1,i]
        m,b=linfits[2,i],linfits[3,i]
        y1,y2=x1*m+b,x2*m+b
        plt.plot([x1,x2],[y1,y2])
    plt.draw()
    
    plt.savefig("test_{0}_{1}.png".format(i,j))

def piecewise_linear_regression(x,y,n=6):
    '''Compute piecewise linear regression on x,y in n equal subsections'''
    output=np.empty((6,n))
    nx=x.size
    step=nx/float(n)
    maxval=-9999
    minval=9999
    for i in range(n):
        lower=i*step
        upper=(i+1)*step-1
        # if we aren't at the top or bottom, use 1/5 more of the step range 
        # to make the regression more stable and consistent between steps 
        # this decreases jumps between adjacent lines/steps
        if i!=0:
            lower-=step/3.0
        if i<n-1:
            upper+=step/3.0
        # compute the regression for the current subset using the scipy stats package
        slope, intercept, r, p, err = stats.linregress(x[lower:upper],y[lower:upper])
        # setup the output
        if i<n-1:
            minval=y[0]
            output[:,i]=[x[i*step],x[(i+1)*step],slope,intercept,minval,maxval]
        else:
            # if we are in the last step, subtract 1 from the step range to avoid out of bounds error?
            maxval=y[-1]
            output[:,i]=[x[i*step],x[-1],slope,intercept,minval,maxval]
    return output

def develop_async(hires=None,lowres=None,isPrecip=True,verbose=True,even_xy=False):
    '''Develop an aynchronous regression between two datasets
    
    optional inputs:
        hires = high resolution input dataset to match = array(ntime x ny x nx)
            default:read_wrf()
        lowres  = low resolution input dataset to be converted = array(ntime x ny x nx)
            default:read_narr()
        isPrecip = boolean, is this variable precip or not (if it is perform regression on the log)
            default:True
        verbose = boolean, print out progress updates
            default:True
    Method
        Calculates the asychronous regression parameters between two
        gridded datasets.  
        For each x,y: 
            sort all points in time by their magnitude. Divide into ten sections
            regress sorted values on each other for each section
    return a structure with original xdata,ydata and regression coefficients
        piecewise linear regression coefficients (plr) are a list of 
            (array(4x10),y_position,x_position)
        (startvalue, endvalue, slope,intercept) for each of the ten segments
        '''
    if hires==None:
        print("Reading WRF data")
        wrf_data=read_wrf()
    else:
        wrf_data=hires
    if lowres==None:
        print("Reading NARR data")
        narr_data=read_narr()
    else:
        narr_data=lowres
    allplrs=list()
    if verbose:
        print("  Generating asynchronous regression")
    # book keeping for a nice progress update
    progress=10
    for i in range(wrf_data.shape[1]):
        # progress update
        if verbose:
            if (float(i)/wrf_data.shape[1])>=(progress/100.0):
                print(str(progress)+"% ",end="")
                sys.stdout.flush()
                progress+=10
        for j in range(wrf_data.shape[2]):
            try:
                if wrf_data[:,i,j].mask[0]:
                    print("has a mask")
                else:
                    print("we shouldn't get here")
            except:
                # Only process this data if it is not masked, assume all times are masked if the first one is
                # sorts, computes log and spaces data "evenly" if requested
                if even_xy:
                    (x,y)=make_xy(narr_data[:,i,j].copy(),wrf_data[:,i,j].copy(),precip=isPrecip,n=500)
                else:
                    (x,y)=(narr_data[:,i,j].copy(),wrf_data[:,i,j].copy())
                    # sort data
                    x.sort()
                    y.sort()
                    if isPrecip:
                        nvals=min(len(np.where(x>0)[0]),len(np.where(y>0)[0]))
                        # subset to only values>0
                        x=x[-nvals:]
                        y=y[-nvals:]
                        # force a minimum value so the log can't explode
                        x[x<1E-5]=1E-5
                        y[y<1E-5]=1E-5
                        # convert to log space
                        x=np.log(x)
                        y=np.log(y)
                # compute the piecewise linear regression
                if len(x)>20:
                    plr=piecewise_linear_regression(x,y)
                    # append results to the list of plrs
                    allplrs.append((plr,i,j))
                    # make plots showing the asynchronous regression at EVERY point in space
                    # test_plot(x,y,plr,i,j); print("{0}, {1}".format(i,j)); sys.stdout.flush()
                else:
                    if len(x)>0:
                        print("\n n={0} : [y={1},x={2}], Obs_max={3} Force_max={4}".format(len(x),i,j,y.max(),x.max()))
                    else:
                        print("\n n={0} : [y={1},x={2}], Obs_all={3} Force_all={4}".format(len(x),i,j,y,x))
    if verbose:
        print("100%")
        sys.stdout.flush()
    return allplrs
    # return Bunch(xdata=narr_data,ydata=wrf_data,plr=allplrs)

def correct_the_top(data,maxval,scale=None):
    tmp=np.where(data>maxval)
    if len(tmp[0])==0:
        return
    if data.max()>500:
        print(data.max())
    if scale==None:scale=0.2*maxval
    errs=data[tmp]-maxval
    data[tmp]=scale*(1-1/(1+errs/scale))+maxval

def correct_the_bottom(data,minval,scale=None):
    tmp=np.where(data<minval)
    if len(tmp[0])==0:
        return
    if data.min()<-100:
        print(data.min())
    
    if scale==None:scale=0.2*minval
    errs=minval-data[tmp]
    data[tmp]=minval-scale*(1-1/(1+errs/scale))

def fill_missing(data,missing=-9999):
    tmp=np.where(data==missing)[0]
    if len(tmp)==0:
        return
    if data.max()<=-999:return
    print("Filling "+str(len(tmp))+" missing values")
    goodvals=np.where(data>-200)[0]
    for i in tmp:
        if i==0:
            curval=data[goodvals[0]]
            data[i]=curval
        elif i==data.shape[0]-1:
            curval=data[goodvals[-1]]
            data[i]=curval
        else:
            nextval=np.where(goodvals>i)[0]
            prevval=np.where(goodvals<i)[0]
            if len(nextval)>0:
                nextval=nextval[0]
            else:nextval=-9999
            if len(prevval)>0:
                prevval=prevval[-1]
            else:prevval=-9999
            curval=-9999
            if prevval==-9999:
                if nextval==-9999:curval=-9999
                else:curval=goodvals[nextval]
            if nextval==-9999:
                if prevval==-9999:curval=-9999
                else:curval=goodvals[prevval]
            #if they are both -9999 then we can't do anything, just assign -9999
            #this shouldnt happen because we return early if there are no values >-999
            if curval==-9999:
                curinterp=(i-prevval)/float(nextval-prevval)
                curval = (curinterp*goodvals[nextval]) + ((1-curinterp)*goodvals[prevval])
            data[i]=curval

def apply_plr(x,plr,vmax=1E3,vmin=-40,precip=True):
    """apply piecewise linear regression to x"""
    n=plr.shape[1]
    output=np.zeros(x.shape,dtype=x.dtype)-9999
    for i in range(n):
        # find the current range of data to work with
        if i==0:
            cur=np.where(x<plr[1,i])
        elif i==n-1:
            cur=np.where(x>=plr[0,i])
        else:
            cur=np.where((x>=plr[0,i])&(x<plr[1,i]))
        # apply the regression slope and offset values
        output[cur]=x[cur]*plr[2,i]+plr[3,i]
    # if these are precip data, need to reverse the log transform
    if precip:
        output=np.exp(output)
        correct_the_top(output,np.exp(vmax))
        output[x<-4]=0
        output[output<0.1]=0
    else:
        correct_the_top(output,vmax,scale=10.0)
        fill_missing(output,missing=-9999)
        correct_the_bottom(output,vmin,scale=10.0)
    return output

def apply_async(data,async,vmax=1E10,vmin=-40,isPrecip=True,verbose=True):
    '''Apply asynchronous regression to entire dataset
    
    data = array(n,ny,nx)
    async= list(plr,i,j)
        plr = piewise linear regression coefficients for data[:,i,j]
            = [start_value, end_value, slope, offset]
        i,j = indices into data, e.g. apply plr to data[:,i,j]
    returns array with the same shape as data after applying regression. 
    '''
    output=np.zeros(data.shape,dtype=data.dtype)
    # book keeping for a progress update
    progress=10
    i=0.0
    # loop through all plr entries in the async list
    for plr in async:
        # progress update
        if verbose:
            if (i/len(async))>=(progress/100.0):
                print(str(progress)+"% ",end="")
                sys.stdout.flush()
                progress+=10
            i+=1
        x=data[:,plr[1],plr[2]]
        # if precip then transform to a log scale
        if isPrecip:
            x[x<1E-20]=1E-20
            x=np.log(x)
        vmax=plr[0][-1,-1]
        vmin=plr[0][-2,0]
        # apply piecewise linear regression to data at [:,i,j] and store in output at [:,i,j]
        output[:,plr[1],plr[2]]=apply_plr(x,plr[0],vmax=vmax,vmin=vmin,precip=isPrecip)
    if verbose:
        print("100%")
        sys.stdout.flush()
    # cleanup just in case, could be needed for precip
    tmp=np.where(~np.isfinite(output))
    output[tmp]=0
    print(len(tmp[0]))
    print(output.max(),output.min())
    return output
