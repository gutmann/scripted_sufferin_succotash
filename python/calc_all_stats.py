#!/usr/bin/env python
import argparse
import traceback
import sys
import os
import re
import glob

import numpy as np

from stat_down import stats,driver,agg
import mygis as io
import date_fun
from nc_reader import NC_Reader
from bunch import Bunch

from flushprint import Flushfile
sys.stdout=Flushfile(sys.stdout)


def load_huc_data(hucfilename,geosubset):
    """Load a grid of HUC data from hucfilename within geosubset bounds"""
    
    # use the NC_Reader class to facilitate subseting from lat/lon data
    hucreader=NC_Reader(None,filelist=[hucfilename],
                    geo_subset=geosubset,
                    readvars=["data"],ntimes=1,
                    latvar="lat",lonvar="lon")
    # just read the first element from the reader
    hucdata=hucreader.next()[0]
    # close the reader now that we are done with it. 
    hucreader.close()
    return hucdata


def load_data(files,varname,extra,res,minval=-100,maxval=500):
    """Load data for varname in from all files matching the geographic area defined in extra
    
    files = a list of filenames to load
    varname = the name of the variable to load from those files
    extra = an arbitrary list of information
        [2]=geographic subset to use  (e.g. [latmin,latmax,lonmin,lonmax])
            geosubset=[25,53,-124.7,-67] #CONUS
            geosubset=[35,43,-112.8,-101.7] # subdomain
        [4]=hucfile to aggregate to (e.g. huc2 data on a grid)
        [5]=geographic file to match lat/lon from (e.g. WRF grid vs Maurer grid)
    
    Returns a structure containing the gridded dataset and the HUC aggregated data if appropriate    
    """
    
    hucfilename=extra[4]
    geosubset=extra[2]
    geo_matchfile=extra[5]
    files.sort()
    # if we don't have a geomatchfile then use nearest neighbor resampling
    usenn=(geo_matchfile==None)
    # use the NC_Reader class to facilitate grid matching and subsetting as necessary
    reader=NC_Reader(None,filelist=files,
                    nn=usenn, bilin=(not usenn),
                    geo_subset=geosubset,
                    geomatch_file=geo_matchfile,
                    readvars=[varname],ntimes=365,
                    latvar="lat",lonvar="lon",
                    glatvar="lat",glonvar="lon")
                    
    # load huc data before actually reading the data in
    outputdata=[]
    hucnames=[]
    if hucfilename:
        # if there are more than one hucfile to match iterate over them
        if type(hucfilename==list):
            hucdata=[]
            for f in hucfilename:
                if f[-3:]!=".nc":
                    f+=res+".nc"
                hucdata.append(load_huc_data(f,geosubset))
                hucnames.append(f.split('/')[-1].split("_")[0])
        # else just load this hucfile
        else:
            hucdata=[load_huc_data(hucfilename,geosubset)]
            f=hucfilename
            if f[-3:]!=".nc":
                f+=res+".nc"
            hucnames.append(f.split('/')[-1].split("_")[0])
    
    # now read the data in from the reader
    for data in reader:
        shp=data[0].shape
        data=np.ma.array(data,mask=( (~np.isfinite(data)) | (data>1E10)))
        outputdata.append(data[0].reshape((1,shp[0],shp[1])))
    outputdata=np.concatenate(outputdata,axis=0)
    reader.close()
    
    # remove any possible bad values (e.g. -9999)
    cleanup(outputdata,minval=minval,maxval=maxval)

    # finally aggreate to HUC shapes if necessary
    aggoutdata=[]
    huclist=[]
    if hucfilename:
        for h in hucdata:
            print("aggregating...")
            if varname=="pr":
                minval=0
            else:
                minval=-100
            aggdata=agg.data2hucs(outputdata,h,minvalue=minval,maxvalue=500)
            curdata=np.dstack(aggdata.data)
            aggoutdata.append(curdata.reshape((curdata.shape[1],1,-1)))
            print(aggoutdata[-1].shape)
            huclist.append(aggdata.hucs)
        
    print(outputdata.shape)
    # return a structure containing the data on the original grid (data) 
    # and the huc shapes (aggdata) as well as the corresponding list of hucs
    return Bunch(data=outputdata,aggdata=aggoutdata,hucs=huclist,hucnames=hucnames)

def calc_year_start_dates(filenames):
    """Calculate the starting points of each year for all files in filenames
    
    based on the filename, calculate the year of the first file and go from there. 
    A better way to do this would be to read the actual time data from the files (if/when it exists)
    or at least read the number of data points in each file, but that is tricky with monthly vs annual files. """
    
    # just a stupid catch for WRF data that start in Oct. 1 instead of Jan 1
    if re.match("added_regridded_daily_.*",filenames[0].split("/")[-1]):
        year1=np.float(filenames[0].split(".")[-2][-4:])
        wrfoffset=92
        return np.floor(np.arange((year1%4.0)/4,365*len(filenames),365.25))+wrfoffset
    try:
        year1=filenames[0].split(".")[-2]
        if len(year1)>4:
            year1=year1.split("_")[-2]
        year1=float(year1)-1
    except IndexError:
        year1=float(filenames[0].split(".")[-3])-1
    except ValueError:
        year1=filenames[0].split(".")[-2]
        if len(year1)>4:
            year1=year1.split("_")[-1]
        year1=float(year1)-1
        
    
    return np.floor(np.arange((year1%4.0)/4,365*len(filenames),365.25))

def temp_stats(names,data1,data2,info):
    """Calculate temperature statistics
    
    names = a list of dataset names
    data1 = a list of tmax datasets
    data2 = a list of tmin datasets
    info = a structure with 
            output_base = output filename 'prefix' 
                if it includes 'annual' then additional annual statistics are calculated
            year_starts = indicies into tmin/tmax for the starting point of each year
                    used to calculate e.g. internanual variability, growing season length, etc"""
    for n,tmax,tmin in zip(names,data1,data2):
        out=info.output_base+"_"+n
        tave=(tmax+tmin)/2
        dtr=tmax-tmin
        print(out)
        
        for thisvar,t in zip(["tmin","tmax","tave","dtr"],[tmin,tmax,tave,dtr]):
            print("MAT "+thisvar)
            mean_annual_temp=stats.mean(t,nyears=t.shape[0])
            print(mean_annual_temp.mean())
            io.write(out+"_mean_"+thisvar,mean_annual_temp)
    
            print("interannual "+thisvar)
            interannual=stats.interannual(t,info.year_starts,fun=np.mean)
            print(interannual.mean())
            io.write(out+"_interannual_"+thisvar,interannual)

            hist=np.vstack(stats.histogram(t,precip=False))
            io.write(out+"_histogram_"+thisvar,hist)
        
        if re.match(".*annual",n):
            frostdays,growing_season=stats.temperature_indicies(tmin,info.year_starts)
            print("Growing season")
            print(growing_season.mean())
            io.write(out+"_growing_season",growing_season)
            print("frostdays")
            print(frostdays.mean())
            io.write(out+"_frostdays",frostdays)

def precip_stats(names,data,info):
    """Calculate precipitation statistics
    
    names = a list of dataset names
    data = a list of precipitation datasets
    info = a structure with 
            output_base = output filename 'prefix' 
                if it includes 'annual' then additional annual statistics are calculated
            year_starts = indicies into tmin/tmax for the starting point of each year
                    used to calculate e.g. internanual variability, growing season length, etc
            precip_threshold = minimum threshold to count as a wet day (e.g. 0, 1, 0.1mm)
            period = length of time period being evaluated in days e.g. for annual=365.25 for monthly~=30
         for extreme events:
            nday_lengths = list of days in aggregation periods for extreme events e.g. 1,2,3 day total
            distributionname = a string for the type of distribution function to use calculating extreme event
                'gamma', 'exponential', or 'weibull' """
    for n,d in zip(names,data):
        out=info.output_base+"_"+n
        print(out)
        wetdaystats=stats.wetdry(d,threshold=info.precip_threshold)
        print("wetfrac,wetspell,dryspell")
        for w in wetdaystats:print(w.mean())
        io.write(out+"_wetfrac",wetdaystats[0])
        io.write(out+"_wetspell",wetdaystats[1])
        io.write(out+"_dryspell",wetdaystats[2])
    
        print("MAP")
        mean_annual_precip=stats.mean(d,nyears=d.shape[0]/info.period)
        print(mean_annual_precip.mean())
        io.write(out+"_MAP",mean_annual_precip)
    
        print("interannual")
        interannual=stats.interannual(d,info.year_starts,fun=np.mean)*info.period
        print(interannual.mean())
        io.write(out+"_interannual",interannual)

        hist=np.vstack(stats.histogram(d))
        io.write(out+"_histogram",hist)
        
        print("WARNING: NOT CALCULATING EXTREME EVENTS")
        print("WARNING: NOT CALCULATING EXTREME EVENTS")
        # commented out for now because extreme events take a long time to calculate
        # if re.match(".*annual",n):
        #     for ndays in info.nday_lengths:
        #         data2test=d[ndays:,...].copy()
        #         for i in range(ndays):
        #             data2test+=d[i:-(ndays-i),...]
        #         print("extremes "+str(ndays+1))
        #         extremes=stats.extremes(data2test,dist_name=info.distributionname)
        #         if extremes!=None:
        #             for i in range(extremes.shape[0]):
        #                 cur=extremes[i,...]
        #                 print(cur[cur>0].mean())
        #             io.write(out+"_extremes_nday"+str(ndays+1),extremes)

def cleanup(data,minval=-999,maxval=1e5):
    sz=data.shape
    for i in range(1,sz[0]):
        tmp=np.where((data[i,...]<minval)|(data[i,...]>maxval)|(~np.isfinite(data[i,...])))
        if len(tmp[0]):
            data[i,tmp[0],tmp[1]]=data[i-1,tmp[0],tmp[1]]
    
def calc_dates(files,ntimes):
    """Calcualte the corresponding dates for a given dataset
    
    Kind of a hack around multiple file / time / date structures
    """
    # for WRF data
    if re.match("added_regridded_daily_.*",files[0].split("/")[-1]):
        year1=np.float(files[0].split(".")[-2][-4:])
        mjd1=date_fun.date2mjd(year1-1,10,1,12,0,0)
        mjd=mjd1+np.arange(ntimes)
        dates=date_fun.mjd2date(mjd)
        return Bunch(year=dates[:,0],month=dates[:,1],day=dates[:,2])
    # for one set of SD data
    try:
        year1=files[0].split(".")[-2]
        if len(year1)>4:
            year1=year1.split("_")[-2]
        year1=float(year1)-1
    # for another set of SD data
    except IndexError:
        year1=files[0].split(".")[-3]
    except ValueError:
        year1=files[0].split(".")[-2]
        if len(year1)>4:
            year1=year1.split("_")[-1]
        year1=float(year1)-1
    
    # once we have the starting year, calculate all other years as modified julian day
    # and convert back to Year, Month, Day dates
    year1=int(year1)
    mjd1=date_fun.date2mjd(year1,1,1,12,0,0)
    mjd=mjd1+np.arange(ntimes)
    dates=date_fun.mjd2date(mjd)
    return Bunch(year=dates[:,0],month=dates[:,1],day=dates[:,2])

def calc_stats(files,v,output_base,info,extra):
    """Calculate primary statistics for a given meteorologic dataset"""
    
    # parse the extra list into a meaningful datastructure
    metadata=Bunch(output_base=extra[3]+output_base,
                   resolution=info[3],
                   distributionname=extra[0],
                   precip_threshold=extra[1],
                   nday_lengths=np.arange(5),
                   year_starts=calc_year_start_dates(files),
                   period=365.25)

    print("Loading Data")
    # load data from files *assumes you can store all data in memory*
    # specifies different valid ranges for different variables
    if v=="pr":
        data=load_data(files,v,extra,metadata.resolution,minval=0,maxval=300)
    else:
        data=load_data(files,v,extra,metadata.resolution,minval=-60,maxval=100)
    # convert data structures to lists so that HUC and full_res can be processed together
    alldata=data.aggdata
    alldata.append(data.data)
    allnames=data.hucnames
    allnames.append("full_res")
    if v=="tasmax":
        # if we just loaded tmax we also need to load tmin
        # first convert filenames
        files2=[f.replace("tasmax","tasmin") for f in files]
        files2=[f.replace("tmax","tmin") for f in files2]
        # now load data2 (tmin)
        data2=load_data(files2,"tasmin",extra,metadata.resolution,minval=-60,maxval=100)
        # similarly aggregate full_res and HUC datasets
        alldata2=data2.aggdata
        alldata2.append(data2.data)
        
        # enforce Tmax>Tmin if not reverse them
        # this may actually make data look better than they are, but it is common practice. 
        for a1,a2 in zip(alldata,alldata2):
            tmp=np.where(a1<a2)
            if len(tmp[0])>0:
                bada1=a1[tmp].copy()
                a1[tmp]=a2[tmp]
                a2[tmp]=bada1
        
    dates=calc_dates(files,alldata[0].shape[0])
    
    # ---------- Calculate Annual values --------------
    fullyearnames=[oldname+"_annual" for oldname in allnames]
    if v=="pr":
        precip_stats(fullyearnames,alldata,metadata)
    else:
        temp_stats(fullyearnames,alldata,alldata2,metadata)
    
    # ---------- Calculate Monthly values --------------
    monthlength=[31,28.25,31,30,31,30,31,31,30,31,30,31]
    for month in range(12):
        # subset data to the current month we are processing
        curtimes=np.where(dates.month==month+1)[0]
        newnames=[oldname+"_month{0:02}".format(month+1) for oldname in allnames]
        newdata=[olddata[curtimes,...] for olddata in alldata]
        metadata.period=monthlength[month]
        metadata.year_starts=[0]
        metadata.year_starts.extend(list(np.where(np.diff(curtimes)>1)[0]+1))
        print(metadata.year_starts)
        
        if v=="pr":
            precip_stats(newnames,newdata,metadata)
        else:
            newdata2=[olddata[curtimes,...] for olddata in alldata2]
            temp_stats(newnames,newdata,newdata2,metadata)

    # ---------- Calculate Seasonal values --------------
    dates.month=(dates.month+1)%12
    for season in range(1,12,3):
        curtimes=np.where(np.abs(dates.month-season)<=1)[0]
        newnames=[oldname+"_season{0}".format(season/3+1) for oldname in allnames]
        newdata=[olddata[curtimes,...] for olddata in alldata]
        metadata.period=monthlength[(season+10)%12]+monthlength[season-1]+monthlength[season]
        metadata.year_starts=[0]
        metadata.year_starts.extend(list(np.where(np.diff(curtimes)>1)[0]+1))
        print(metadata.year_starts)
        
        if v=="pr":
            precip_stats(newnames,newdata,metadata)
        else:
            newdata2=[olddata[curtimes,...] for olddata in alldata2]
            temp_stats(newnames,newdata,newdata2,metadata)

def mk_all_dirs(dirname):
    """Loop over all subdirectories in dirname and create any missing directories"""
    original_dir=os.getcwd()
    subdirs=dirname.split("/")
    for curdir in subdirs:
        if os.path.isfile(curdir):
            os.rename(curdir,os.tempnam("./",curdir+"_"))
        if not os.path.isdir(curdir):
            try:
                os.mkdir(curdir)
            except Exception as e:
                print(e)
        os.chdir(curdir)

    os.chdir(original_dir)


if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Calculate nday/nyear return period events. ')
        parser.add_argument('-method',dest="methods",nargs="?",action='store',help="SD methods to test [SDmon,CAe0,SARe0,SDe0,CA]",
                    default=["SDmon","CAe0","SARe0","SDe0","CA","CAe1","SDe1"])#,"SARe1","SDe","CAe","SD"])
        parser.add_argument('-model',dest="models",nargs="?",action='store',
                    default=["ncep","narr"],help="model forcing to test [ccsm,ncep,narr]")
        parser.add_argument('-res',dest="resolution",nargs="?",action='store',help="resolution to run [12km,6km]",
                    default=["12km","6km"])
        parser.add_argument('-var',dest="variable",nargs="?",action='store',
                    default=["pr"], help="variable to test ([pr],tasmax)")
        parser.add_argument('-bc',dest="BC",nargs="?",action='store',help="Bias Corrected, or not [BC],nobc ",
                    default=["BC"])#,""])
        # bd="/glade/scratch/gutmann/usbr/hucdata/"
        bd="/d2/gutmann/usbr/hucdata/"
        parser.add_argument('-huc',dest="huc",nargs="?",action='store',help="HUC file names [default= all]",
                    default=[bd+"HUC02/huc2_",bd+"HUC04/huc4_",bd+"HUC08/huc8_"])#,bd+"HUC12/huc12_"])
        parser.add_argument('-out',dest="outputdir",nargs="?",action='store',
                    default="./",help="output directory [./]")
        parser.add_argument('-yearsearch',dest="yearsearch",nargs="?",action='store',
                    default=["200[1-8]"],help='years to search for [200[1-8]]')
                    # default=["19*"],help='years to search for [200[1-8]]')
        parser.add_argument ('-distribution', nargs="?", action='store',
                default="gamma", help='name of extreme distribution ([gamma],weibull,exponential)', dest='distribution')
        parser.add_argument ('-pr_threshold', nargs="?", action='store', dest='prthresh',
                default="0.0", help='Threshold to use when calculating wet day status')
        parser.add_argument ('--runsub', action='store_true',
                default=False, help='Run subdomain only [False]', dest='runsub')
        parser.add_argument ('--runobs', action='store_true',
                default=False, help='Run observations only [False]', dest='runobs')
        parser.add_argument ('--runwrf', action='store_true',
                default=False, help='Run WRF only [False]', dest='runwrf')
        parser.add_argument ('--runforcing', action='store_true',
                default=False, help='Run Forcing only [False]', dest='runforcing')
        
        parser.add_argument('-v', '--version',action='version',
                version='calc_all_stats v.1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        distribution=args.distribution
        outputdir=args.outputdir
        hucfile=args.huc
        # hucfile=None
        pr_threshold=float(args.prthresh)
        print("precip threshold="+str(pr_threshold))
        
        runobs=args.runobs
        runforcing=args.runforcing
        runwrf=args.runwrf
        runstat=(not runobs) and (not runforcing) and (not runwrf)
        
        if args.BC=="nobc":args.BC=""
        geosubset=[25,53,-124.7,-67] #CONUS
        if args.runsub:
            geosubset=[35,43,-112.8,-101.7] # subdomain
        
        # driver.drive requires lists to iterate over, but CLI args will be individual
        #  elements.  However, default values are all lists, so we don't want to make them
        # lists of lists so we have to test to see if they are a list already first...
        for k in args.__dict__.keys():
            if type(args.__dict__[k])!=list:
                args.__dict__[k]=[args.__dict__[k]]
        if outputdir[-1]!='/':outputdir+='/'
        if outputdir!="./":
            mk_all_dirs(outputdir[:-1])
            
        georeffile=None
        print(geosubset)
        if runobs:
            print("Running Observations Only")
            sys.stdout.flush()
        if runforcing:
            hucfile=None
            print("Running Forcing Only")
            georeffile="../obs/maurer.125/pr/nldas_met_update.obs.daily.pr.2000.nc"
            sys.stdout.flush()
        
        # 
        if runwrf:
            geosubset=[35,43,-112.8,-101.7]
            files=glob.glob("../wrf/4km/added*")
            files.sort()
            calc_stats(files,"pr","wrf",[0,0,0,"12km"],["gamma",pr_threshold,geosubset,outputdir,hucfile,None])
            os._exit(0)
        
        
        # driver iterates over all combinations of methods/variables/forcing/resolutions and calls calc_stats for each
        driver.drive(calc_stats, #function to call
                    yearsearch=args.yearsearch[0], #subset of years to run specified as a glob expression
                    obs=runobs,stat=runstat,runforce=runforcing, # options to run either SD methods, observations or forcing
                    extra=[distribution,pr_threshold,geosubset,outputdir,hucfile,georeffile], #passed to calc_stats
                    methods=args.methods,   # list of SD methods to run
                    BCs=args.BC,            # list of BC/noBC
                    models=args.models,     # list of forcing models ([NCEP,NARR])
                    resolutions=args.resolution,# list of resolutions ([6km,12km])
                    variables=args.variable,    # list of variables ([pr,tasmax,tasmin])
                    extendedmethods=False)  #extendedmethods includes hot,dry,wet,cold experiments (don't use)
                    
    # overarching error handling to print out something useful if it crashes
    except Exception as e:
        print("Error Unexpected Exception")
        print(e)
        traceback.print_exc()
        os._exit(1)
