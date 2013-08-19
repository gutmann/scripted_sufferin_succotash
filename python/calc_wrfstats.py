#!/usr/bin/env python

# Note, this file is only minor tweaks from calc_all_stats.py specific to WRF data
# for questions read calc_all_stats

import argparse
import traceback
import sys
import os
import re
import glob

import numpy as np

from stat_down import stats,driver,agg
from stat_down import myio as io
import date_fun
from nc_reader import NC_Reader
from bunch import Bunch

def load_huc_data(hucfilename,geosubset):
    hucreader=NC_Reader(None,filelist=[hucfilename],
                    geo_subset=geosubset,
                    readvars=["data"],ntimes=1,
                    latvar="lat",lonvar="lon")
    hucdata=hucreader.next()[0]
    hucreader.close()
    return hucdata


def load_data(files,varname,extra,res,minval=-100,maxval=500):
    
    hucfilename=extra[4]
    geosubset=extra[2]
    geo_matchfile=extra[5]
    files.sort()
    usenn=(geo_matchfile==None)
    reader=NC_Reader(None,filelist=files,
                    nn=usenn, bilin=(not usenn),
                    geo_subset=geosubset,
                    geomatch_file=geo_matchfile,
                    readvars=[varname],ntimes=365,
                    latvar="lat",lonvar="lon",
                    glatvar="lat",glonvar="lon")
    outputdata=[]
    hucnames=[]
    hucfilename=None
    if hucfilename:
        if type(hucfilename==list):
            hucdata=[]
            for f in hucfilename:
                if f[-3:]!=".nc":
                    f+=res+".nc"
                hucdata.append(load_huc_data(f,geosubset))
                hucnames.append(f.split('/')[-1].split("_")[0])
        else:
            hucdata=[load_huc_data(hucfilename,geosubset)]
            f=hucfilename
            if f[-3:]!=".nc":
                f+=res+".nc"
            hucnames.append(f.split('/')[-1].split("_")[0])
        
    huclist=[]
    for data in reader:
        shp=data[0].shape
        outputdata.append(data[0].reshape((1,shp[0],shp[1])))
    outputdata=np.concatenate(outputdata,axis=0)
    reader.close()
    # remove any possible bad values (e.g. -9999)
    cleanup(outputdata,minval=minval,maxval=maxval)

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
    return Bunch(data=outputdata,aggdata=aggoutdata,hucs=huclist,hucnames=hucnames)

def calc_year_start_dates(filenames):
    if re.match("added_regridded_daily_.*",filenames[0].split("/")[-1]):
        year1=np.float(filenames[0].split(".")[-2][-4:])
        wrfoffset=92
        return np.floor(np.arange((year1%4.0)/4,365*len(filenames),365.25))+wrfoffset
    if re.match("daily_.*",filenames[0].split("/")[-1]):
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
    
    return np.floor(np.arange((year1%4.0)/4,365*len(filenames),365.25))

def temp_stats(names,data1,data2,info):
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
        
        if re.match(".*annual",n):
            for ndays in info.nday_lengths:
                data2test=d[ndays:,...].copy()
                for i in range(ndays):
                    data2test+=d[i:-(ndays-i),...]
                print("extremes "+str(ndays+1))
                extremes=stats.extremes(data2test,dist_name=info.distributionname)
                if extremes!=None:
                    for i in range(extremes.shape[0]):
                        cur=extremes[i,...]
                        print(cur[cur>0].mean())
                    io.write(out+"_extremes_nday"+str(ndays+1),extremes)

def cleanup(data,minval=-999,maxval=1e5):
    sz=data.shape
    for i in range(1,sz[0]):
        tmp=np.where((data[i,...]<minval)|(data[i,...]>maxval)|(~np.isfinite(data[i,...])))
        if len(tmp[0]):
            data[i,tmp[0],tmp[1]]=data[i-1,tmp[0],tmp[1]]
    
def calc_dates(files,ntimes):
    if re.match("added_regridded_daily_.*",files[0].split("/")[-1]):
        year1=np.float(files[0].split(".")[-2][-4:])
        mjd1=date_fun.date2mjd(year1-1,10,1,12,0,0)
        mjd=mjd1+np.arange(ntimes)
        dates=date_fun.mjd2date(mjd)
        return Bunch(year=dates[:,0],month=dates[:,1],day=dates[:,2])
    if re.match("daily_.*",files[0].split("/")[-1]):
        year1=np.float(files[0].split(".")[-2][-4:])
        mjd1=date_fun.date2mjd(year1-1,10,1,12,0,0)
        mjd=mjd1+np.arange(ntimes)
        dates=date_fun.mjd2date(mjd)
        return Bunch(year=dates[:,0],month=dates[:,1],day=dates[:,2])
    try:
        year1=files[0].split(".")[-2]
        if len(year1)>4:
            year1=year1.split("_")[-2]
        year1=float(year1)-1
    except IndexError:
        year1=files[0].split(".")[-3]
    year1=int(year1)
    mjd1=date_fun.date2mjd(year1,1,1,12,0,0)
    mjd=mjd1+np.arange(ntimes)
    dates=date_fun.mjd2date(mjd)
    return Bunch(year=dates[:,0],month=dates[:,1],day=dates[:,2])

def calc_stats(files,v,output_base,info,extra):
    metadata=Bunch(output_base=extra[3]+output_base,
                   resolution=info[3],
                   distributionname=extra[0],
                   precip_threshold=extra[1],
                   nday_lengths=np.arange(5),
                   year_starts=calc_year_start_dates(files),
                   period=365.25)

    print(metadata.year_starts)
    print("Loading Data")
    if v=="pr":
        data=load_data(files,v,extra,metadata.resolution,minval=0,maxval=300)
    else:
        data=load_data(files,v,extra,metadata.resolution,minval=-60,maxval=100)
    alldata=data.aggdata
    alldata.append(data.data)
    allnames=data.hucnames
    allnames.append("full_res")
    if v=="tasmax":
        files2=[f.replace("tasmax","tasmin") for f in files]
        files2=[f.replace("tmax","tmin") for f in files2]
        data2=load_data(files2,"tasmin",extra,metadata.resolution,minval=-60,maxval=100)
        alldata2=data2.aggdata
        alldata2.append(data2.data)
        # enforce Tmax>Tmin if not reverse them
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
    for i in range(len(subdirs)):
        os.chdir("../")

if __name__ == '__main__':
    
    pr_threshold=0
    geosubset=None
    outputdir="./"
    bd="/glade/scratch/gutmann/usbr/hucdata/"
    hucfile=[bd+"HUC02/huc2_",bd+"HUC04/huc4_",bd+"HUC08/huc8_"]
    hucfiles=None
    files=glob.glob("/glade/u/home/gutmann/scratch/wrfoutput/4km/daily*")
    calc_stats(files,"pr","wrf",[0,0,0,"12km"],["gamma",pr_threshold,geosubset,outputdir,hucfile,None])
    
    # try:
    #     parser= argparse.ArgumentParser(description='Calculate nday/nyear return period events. ')
    #     parser.add_argument('-method',dest="methods",nargs="?",action='store',help="SD methods to test [SDmon,CAe0,SARe0,SDe0,CA]",
    #                 default=["SDmon","CAe0","SARe0","SDe0","CA","CAe1","SDe1"])#,"SARe1","SDe","CAe","SD"])
    #     parser.add_argument('-model',dest="models",nargs="?",action='store',
    #                 default=["ncep","narr"],help="model forcing to test [ncep,narr]")
    #     parser.add_argument('-res',dest="resolution",nargs="?",action='store',help="resolution to run [12km,6km]",
    #                 default=["12km","6km"])
    #     parser.add_argument('-var',dest="variable",nargs="?",action='store',
    #                 default=["pr"], help="variable to test ([pr],tasmax)")
    #     parser.add_argument('-bc',dest="BC",nargs="?",action='store',help="Bias Corrected, or not [BC],nobc ",
    #                 default=["BC"])#,""])
    #     bd="/glade/scratch/gutmann/usbr/hucdata/"
    #     parser.add_argument('-huc',dest="huc",nargs="?",action='store',help="HUC file names [default= all]",
    #                 default=[bd+"HUC02/huc2_",bd+"HUC04/huc4_",bd+"HUC08/huc8_"])#,bd+"HUC12/huc12_"])
    #     parser.add_argument('-out',dest="outputdir",nargs="?",action='store',
    #                 default="./",help="output directory [./]")
    #     parser.add_argument('-yearsearch',dest="yearsearch",nargs="?",action='store',
    #                 default=["200[1-8]"],help='years to search for [200[1-8]]')
    #                 # default=["19*"],help='years to search for [200[1-8]]')
    #     parser.add_argument ('-distribution', nargs="?", action='store',
    #             default="gamma", help='name of extreme distribution ([gamma],weibull,exponential)', dest='distribution')
    #     parser.add_argument ('-pr_threshold', nargs="?", action='store', dest='prthresh',
    #             default="0.0", help='Threshold to use when calculating wet day status')
    #     parser.add_argument ('--runobs', action='store_true',
    #             default=False, help='Run observations only [False]', dest='runobs')
    #     parser.add_argument ('--runwrf', action='store_true',
    #             default=False, help='Run WRF only [False]', dest='runwrf')
    #     parser.add_argument ('--runforcing', action='store_true',
    #             default=False, help='Run Forcing only [False]', dest='runforcing')
    #     
    #     parser.add_argument('-v', '--version',action='version',
    #             version='calc_all_stats v.1.0')
    #     parser.add_argument ('--verbose', action='store_true',
    #             default=False, help='verbose output', dest='verbose')
    #     args = parser.parse_args()
    # 
    #     distribution=args.distribution
    #     outputdir=args.outputdir
    #     hucfile=args.huc
    #     # hucfile=None
    #     pr_threshold=float(args.prthresh)
    #     print("precip threshold="+str(pr_threshold))
    #     
    #     runobs=args.runobs
    #     runforcing=args.runforcing
    #     runwrf=args.runwrf
    #     runstat=(not runobs) and (not runforcing) and (not runwrf)
    #     
    #     if args.BC=="nobc":args.BC=""
    #     geosubset=[25,53,-124.7,-67]
    #     # geosubset=[35,43,-113,-101]
    #     # geosubset=[35,43,-112.8,-101.7]
    #     # driver.drive requires lists to iterate over, but CLI args will be individual
    #     #  elements.  However, default values are all lists, so we don't want to make them
    #     # lists of lists so we have to test to see if they are a list already first...
    #     for k in args.__dict__.keys():
    #         if type(args.__dict__[k])!=list:
    #             args.__dict__[k]=[args.__dict__[k]]
    #     if outputdir[-1]!='/':outputdir+='/'
    #     if outputdir!="./":
    #         mk_all_dirs(outputdir[:-1])
    #         
    #     georeffile=None
    #     print(geosubset)
    #     if runobs:
    #         print("Running Observations Only")
    #         sys.stdout.flush()
    #     if runforcing:
    #         hucfile=None
    #         print("Running Forcing Only")
    #         georeffile="../obs/maurer.125/pr/nldas_met_update.obs.daily.pr.2000.nc"
    #         sys.stdout.flush()
    #         
    #     if runwrf:
    #         geosubset=[35,43,-112.8,-101.7]
    #         files=glob.glob("/glade/u/home/gutmann/scratch/wrfoutput/4km/added*")
    #         files.sort()
    #         calc_stats(files,"pr","wrf",[0,0,0,"12km"],["gamma",pr_threshold,geosubset,outputdir,hucfile,None])
    #         os._exit(0)
    #         
    #     driver.drive(calc_stats,
    #                 yearsearch=args.yearsearch[0],
    #                 obs=runobs,stat=runstat,runforce=runforcing,
    #                 extra=[distribution,pr_threshold,geosubset,outputdir,hucfile,georeffile],
    #                 methods=args.methods,
    #                 BCs=args.BC,
    #                 models=args.models,
    #                 resolutions=args.resolution,
    #                 variables=args.variable,
    #                 extendedmethods=False)
    #                 
    # except Exception as e:
    #     print("Error Unexpected Exception")
    #     print(e)
    #     traceback.print_exc()
    #     os._exit(1)
