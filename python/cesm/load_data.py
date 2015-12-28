import numpy as np
import mygis
import glob

from bunch import Bunch

current_precc_file="b.e11.B20TRC5CNBDRD.f09_g16.{:03}.cam.h0.PRECC.192001-200512.nc"
current_precl_file="b.e11.B20TRC5CNBDRD.f09_g16.{:03}.cam.h0.PRECL.192001-200512.nc"
future_precc_file="b.e11.BRCP85C5CNBDRD.f09_g16.{:03}.cam.h0.PRECC.200601-*.nc"
future_precl_file="b.e11.BRCP85C5CNBDRD.f09_g16.{:03}.cam.h0.PRECL.200601-*.nc"

current_start_year=1920
future_start_year=2006


y0=130
y1=145
x0=200
x1=208

days_per_month=[31,28,31,30,31,30,31,31,30,31,30,31]

def monthly_means(verbose=False):
    """docstring for main"""
    current_start_point = (1990-current_start_year)*12
    current_end_point = (2000-current_start_year)*12

    future_start_point = (2025-future_start_year)*12
    future_end_point = (2035-future_start_year)*12
    months=np.arange(1,13)
    
    cur=[]
    fut=[]
    delta=[]
    
    for i in range(30):
        if verbose:print(i)
        filenum=i+1
        if filenum>35:filenum+=100-35
        currentdata=(mygis.read_nc(current_precc_file.format(filenum),"PRECC",
                    returnNCvar=True).data[current_start_point:current_end_point, y0:y1,x0:x1].mean(axis=1).mean(axis=1)
                    +mygis.read_nc(current_precl_file.format(filenum),"PRECL",
                    returnNCvar=True).data[current_start_point:current_end_point, y0:y1,x0:x1].mean(axis=1).mean(axis=1))
        futuredata =(mygis.read_nc(future_precc_file.format(filenum), "PRECC",
                    returnNCvar=True).data[future_start_point:future_end_point,   y0:y1,x0:x1].mean(axis=1).mean(axis=1)
                    +mygis.read_nc(future_precl_file.format(filenum), "PRECL",
                    returnNCvar=True).data[future_start_point:future_end_point,   y0:y1,x0:x1].mean(axis=1).mean(axis=1))
        
        cur.append(1000.*86400.*np.reshape(currentdata,(10,12)).mean(axis=0))
        fut.append(1000.*86400.*np.reshape(futuredata,(10,12)).mean(axis=0))
        delta.append(fut-cur)
    
    return Bunch(current=cur, future=fut, delta=delta, times=months)
        
        
def monthly_timeseries(verbose=False):
    """docstring for main"""
    
    output=[]
    yearly=[]
    ondjfma_data=[]
    warm_data=[]
    
    ens_mean_yearly=None
    ens_mean=None
    global current_precc_file
    global current_precl_file
    global current_start_point
    global current_end_point
    n_ens=40.0
    
    for i in range(int(n_ens)):
        current_start_point=0
        current_end_point=None
        future_start_point=0
        future_end_point=900
        
        filenum=i+1
        if filenum>35:filenum+=100-35
        if verbose:print("Ensemble={}".format(filenum))
        if not glob.glob(current_precc_file.format(filenum)):
            current_precc_file=current_precc_file.replace("192001-","185001-")
            current_precl_file=current_precl_file.replace("192001-","185001-")
            start_is_1850=True
            current_start_point+=840
            # current_end_point+=840
        else:
            start_is_1850=False
            
        if verbose:
            print(i)
            print(current_precc_file.format(filenum))
            print(future_precc_file.format(filenum))
            
        currentdata=(mygis.read_nc(current_precc_file.format(filenum),"PRECC",
                    returnNCvar=True).data[current_start_point:current_end_point, y0:y1,x0:x1].mean(axis=1).mean(axis=1)
                    +mygis.read_nc(current_precl_file.format(filenum),"PRECL",
                    returnNCvar=True).data[current_start_point:current_end_point, y0:y1,x0:x1].mean(axis=1).mean(axis=1))
        futuredata =(mygis.read_nc(future_precc_file.format(filenum), "PRECC",
                    returnNCvar=True).data[future_start_point:future_end_point,   y0:y1,x0:x1].mean(axis=1).mean(axis=1)
                    +mygis.read_nc(future_precl_file.format(filenum), "PRECL",
                    returnNCvar=True).data[future_start_point:future_end_point,   y0:y1,x0:x1].mean(axis=1).mean(axis=1))
        
        if verbose:
            print(currentdata.shape)
            print(futuredata.shape)
        
        if start_is_1850:
            current_precc_file=current_precc_file.replace("185001-","192001-")
            current_precl_file=current_precl_file.replace("185001-","192001-")
            current_start_point-=840
            
        
        data=np.concatenate([currentdata,futuredata],axis=0)
        for month in range(data.shape[0]):
            data[month]*=1000.0*86400.0*days_per_month[month%12]
        
        nyears = data.shape[0]/12
        annual = data[:nyears*12].reshape((nyears,12)).sum(axis=1)
        ondjfma= np.take(data[:nyears*12].reshape((nyears,12)), [0,1,2,3,9,10,11],axis=1).sum(axis=1)
        warm   = np.take(data[:nyears*12].reshape((nyears,12)), [4,5,6,7,8],axis=1).sum(axis=1)
        
        output.append(data)
        yearly.append(annual)
        ondjfma_data.append(ondjfma)
        warm_data.append(warm)
        if ens_mean==None:
            ens_mean        = np.zeros(data.shape)+data
            ens_mean_yearly = np.zeros(annual.shape)+annual
            ens_mean_ondjfma= np.zeros(ondjfma.shape)+ondjfma
            ens_mean_warm   = np.zeros(warm.shape)+warm
        else:
            n=ens_mean.shape[0]
            ens_mean+=data[:n]
            n=ens_mean_yearly.shape[0]
            ens_mean_yearly+=annual[:n]
            n=ens_mean_ondjfma.shape[0]
            ens_mean_ondjfma+=ondjfma[:n]
            n=ens_mean_warm.shape[0]
            ens_mean_warm+=warm[:n]
            
            
    ens_mean/=n_ens
    ens_mean_yearly/=n_ens
    ens_mean_ondjfma/=n_ens
    ens_mean_warm/=n_ens
    
    times=np.arange(ens_mean.shape[0])/12.0 + current_start_year
    years=np.arange(annual.shape[0]) + current_start_year
    return Bunch(data=output, times=times, mean=ens_mean, yearly=yearly, yrmean=ens_mean_yearly, years=years, 
                 ondjfma=ondjfma_data, coolmean=ens_mean_ondjfma,
                 warm=warm_data, warmmean=ens_mean_warm)
    
        
