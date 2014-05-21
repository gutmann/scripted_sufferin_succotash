from __future__ import print_function
import sys
import numpy as np
from scipy import stats
from bunch import Bunch

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
    if len(x)<24:
        n=int(np.floor(len(x)/4))
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
                        x[x<1E-10]=1E-10
                        y[y<1E-10]=1E-10
                        # convert to log space
                        x=np.log(x)
                        y=np.log(y)
                # compute the piecewise linear regression
                if len(x)>=5: #require a minimum number of data points per line segment...
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
    curmax=data.max()
    if (curmax>500) and (curmax<1e5):
        print(curmax,maxval)
    if scale==None:scale=0.2*maxval
    errs=data[tmp]-maxval
    # rescale high data to asymptotically approach maxval+scale, scale is often 10[K] or 20% max[mm]
    data[tmp]=scale*(1-1/(1+errs/scale))+maxval

def correct_the_bottom(data,minval,scale=None):
    tmp=np.where(data<minval)
    if len(tmp[0])==0:
        return
    curmin=data.min()
    if (curmin<-100) and (curmin>-1e5):
        print(curmin)
    
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

def apply_plr(x,plr,vmax=1E4,vmin=-30,precip=True):
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
        try:
            output[cur]=x[cur]*plr[2,i]+plr[3,i]
        except:
            # we get here when plr[2:3,i]==NaN bevause x ~= -11 ~= 1e-5 mm
            output[cur]=0
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

def apply_async(data,async,vmax=1E4,vmin=-30,isPrecip=True,verbose=True):
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
