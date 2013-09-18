#!/usr/bin/env python
from __future__ import print_function
import argparse,sys,os,traceback
from glob import glob
import datetime

import numpy as np
import Nio

import swim_io
from bunch import Bunch

res=36

if res==4:
    wrf_location="/glade/scratch/evanvyve/wrfoutput/tiny_4km/"
    wrfsearch="200[7-8]/*00"
elif res==36:
    wrf_location="/glade/scratch/gutmann/wrfoutput/36km/"
    wrfsearch="wrfout_d01_200[7-8]*00"

def write_file(filename,inputnc,variable):
    """write data from data structure into file"""
    filename.replace("4km","hires")
    filename.replace("36km","lores")
    NCfile=Nio.open_file(filename,"w",format="nc")
    # create all dimensions we need
    NCfile.create_dimension('Time', None)
    NCfile.create_dimension('DateStrLen', 19)
    NCfile.create_dimension('bottom_top', 20)
    NCfile.create_dimension('south_north', 1)
    NCfile.create_dimension('west_east', 1)
    NCfile.create_dimension('south_north_stag', 1)
    NCfile.create_dimension('west_east_stag', 1)
    
    # setup south north and east west dimensions depending on whether we are working with U or V data
    if variable=="U":
        SNdim="south_north"
        EWdim="west_east_stag"
    else:
        SNdim="south_north_stag"
        EWdim="west_east"
    
    # create the variables we need
    NCfile.create_variable("Times","S1",("Time","DateStrLen"))
    NCfile.create_variable(variable,"f",("Time","bottom_top",SNdim,EWdim))
    NCfile.create_variable("XLAT","f",("Time","south_north","west_east"))
    NCfile.create_variable("XLONG","f",("Time","south_north","west_east"))
    NCfile.create_variable("XLAT_U","f",("Time",SNdim,EWdim))
    NCfile.create_variable("XLONG_U","f",("Time",SNdim,EWdim))
    NCfile.create_variable("XLAT_V","f",("Time",SNdim,EWdim))
    NCfile.create_variable("XLONG_V","f",("Time",SNdim,EWdim))
    
    # store the actual data in the file
    ntimes=inputnc.variable.data.shape[0]
    for i in range(ntimes):
        NCfile.variables["Times"][i,:]=str(inputnc.time[i])
    NCfile.variables[variable][:]=inputnc.variable.data
    NCfile.variables["XLAT"][:]=inputnc.lat.data
    NCfile.variables["XLONG"][:]=inputnc.lon.data
    NCfile.variables["XLAT_U"][:]=inputnc.ulat.data
    NCfile.variables["XLONG_U"][:]=inputnc.ulon.data
    NCfile.variables["XLAT_V"][:]=inputnc.vlat.data
    NCfile.variables["XLONG_V"][:]=inputnc.vlon.data
    
    # copy over all variable attributes to the netcdf file
    for k in inputnc.attributes.keys():
        NCfile.__setattr__(k,inputnc.attributes[k])
    for k in inputnc.variable.attributes.keys():
        NCfile.variables[variable].__setattr__(k,inputnc.variable.attributes[k])
    for k in inputnc.lat.attributes.keys():
        NCfile.variables["XLAT"].__setattr__(k,inputnc.lat.attributes[k])
    for k in inputnc.lon.attributes.keys():
        NCfile.variables["XLONG"].__setattr__(k,inputnc.lon.attributes[k])
    for k in inputnc.ulat.attributes.keys():
        NCfile.variables["XLAT_U"].__setattr__(k,inputnc.ulat.attributes[k])
    for k in inputnc.ulon.attributes.keys():
        NCfile.variables["XLONG_U"].__setattr__(k,inputnc.ulon.attributes[k])
    for k in inputnc.vlat.attributes.keys():
        NCfile.variables["XLAT_V"].__setattr__(k,inputnc.vlat.attributes[k])
    for k in inputnc.vlon.attributes.keys():
        NCfile.variables["XLONG_V"].__setattr__(k,inputnc.vlon.attributes[k])
        
    try:
        NCfile.__setattr__("history","Created by wrf2anen.py \nOn: "+str(datetime.datetime.today())+" by: "+os.environ['USER']+"  \noldhistory{"+NCfile.attributes["history"]+"}")
    except KeyError:
        NCfile.__setattr__("history","Created by wrf2anen.py \nOn: "+str(datetime.datetime.today())+" by: "+os.environ['USER'])
        
    NCfile.close()
    
def setup_attributes(filename,variable):
    """Load in all required attributes from the netcdf file into a data structure"""
    # generic structure of structures... of structures
    # with variable.data prepopulated
    outputdata=Bunch(attributes=Bunch(), variable=Bunch(attributes=Bunch()),
                     lat=Bunch(attributes=Bunch()),   lon=Bunch(attributes=Bunch()),
                     ulat=Bunch(attributes=Bunch()), ulon=Bunch(attributes=Bunch()),
                     vlat=Bunch(attributes=Bunch()), vlon=Bunch(attributes=Bunch()))
    
    ncfile=Nio.open_file(filename,"r",format="nc")
    
    for k in ncfile.attributes.keys():
        outputdata.attributes[k]=ncfile.attributes[k]
    
    for localv,ncv in zip(["variable","lat","lon","ulat","ulon","vlat","vlon"],[variable,"XLAT","XLONG","XLAT_U","XLONG_U","XLAT_V","XLONG_V"]):
        thisvar=ncfile.variables[ncv]
        for k in thisvar.attributes.keys():
            outputdata[localv].attributes[k]=thisvar.attributes[k]

    return outputdata


def main(variable,output_location="./",filestart=0,fileend=None,verbose=False):
    """process WRF data for the AnEn Analog Ensemble"""
    # find input files
    files=glob(wrf_location+wrfsearch)
    # put them in order
    files.sort()
    files=files[filestart:fileend]
    print("Nfiles="+str(len(files)))
    # set up output file location
    outputfile_base=output_location+"wrfout_d01_"+str(res)+"km_"+variable

    # define subdomain
    if res==4:
        x1=125
        x2=150
        y1= 60
        y2=160
        # create variables up front
        ntimes=len(files)*24
    elif res==36:
        x1=13
        x2=28
        y1=9
        y2=23
        ntimes=len(files)*24*130
        
    nz=20
    ny=y2-y1
    nx=x2-x1
    
    print("Primary Memory Requirements="+str(ntimes*26/1024.0*ny/1024.0*nx*4.0/1024.)+"GB")
    curdata=np.zeros((ntimes,nz,ny,nx),dtype=np.float32)
    times=np.empty((ntimes,),dtype="S19")
    times[:]=""
    curlat =np.zeros((1,ny,nx),dtype=np.float32)
    curlon =np.zeros((1,ny,nx),dtype=np.float32)
    curulat=np.zeros((1,ny,nx),dtype=np.float32)
    curulon=np.zeros((1,ny,nx),dtype=np.float32)
    curvlat=np.zeros((1,ny,nx),dtype=np.float32)
    curvlon=np.zeros((1,ny,nx),dtype=np.float32)
    
    # loop over all files reading in data
    i=0
    for f in files:
        if verbose:
            print(f)
            sys.stdout.flush()
        nc_data=Nio.open_file(f,format="nc",mode="r")
        ntimes=nc_data.variables[variable].shape[0]
        # read the main data in
        curdata[i:i+ntimes,...]=nc_data.variables[variable][:,:nz,y1:y2,x1:x2]
        if res==36:
            if f!=files[-1]:
                ntimes-=1
            
        # times/character arrays are awkward to read, there is probably a better way? 
        for t in range(ntimes):
            times[i+t]="".join(nc_data.variables["Times"][t,:])
        # if this is the first time through read in the lat/lon variables too
        # they are constant in time so you don't need to read them every time
        if i==0:
            curlat[0,...]= nc_data.variables["XLAT"][0,y1:y2,x1:x2]
            curlon[0,...]= nc_data.variables["XLONG"][0,y1:y2,x1:x2]
            curulat[0,...]=nc_data.variables["XLAT_U"][0,y1:y2,x1:x2]
            curulon[0,...]=nc_data.variables["XLONG_U"][0,y1:y2,x1:x2]
            curvlat[0,...]=nc_data.variables["XLAT_V"][0,y1:y2,x1:x2]
            curvlon[0,...]=nc_data.variables["XLONG_V"][0,y1:y2,x1:x2]
        i=i+ntimes
        nc_data.close()
    
    curdata=curdata[:i,...]
    times=times[:i]
    # now make the lat lon variables have a time axis (constant in time) 
    # not sure this is necessary, Nio might take of this on the right step... 
    ntimes=curdata.shape[0]
    curlat=curlat.repeat(ntimes,axis=0)
    curlon=curlon.repeat(ntimes,axis=0)
    curulat=curulat.repeat(ntimes,axis=0)
    curulon=curulon.repeat(ntimes,axis=0)
    curvlat=curvlat.repeat(ntimes,axis=0)
    curvlon=curvlon.repeat(ntimes,axis=0)
    
    # read in all the attributes from the first netCDF file (and use them for ALL files...)
    ncdataset=setup_attributes(files[0],variable)
    # store the time data
    ncdataset.time=times
    
    # loop through space writing individual files
    for x in range(x1,x2):
        if verbose:
            print("Progress= {:4}%".format(int(round((x-x1)/float(x2-x1)*100))),end="\r")
            sys.stdout.flush()
        for y in range(y1,y2):
            thisy=y-y1
            thisx=x-x1
            # store the current grid point in the output data structure
            ncdataset.variable.data=curdata[:,:,thisy,thisx].reshape((ntimes,20,1,1))
            ncdataset.ulat.data=curulat[:,thisy,thisx].reshape((ntimes,1,1))
            ncdataset.ulon.data=curulon[:,thisy,thisx].reshape((ntimes,1,1))
            ncdataset.vlat.data=curvlat[:,thisy,thisx].reshape((ntimes,1,1))
            ncdataset.vlon.data=curvlon[:,thisy,thisx].reshape((ntimes,1,1))
            ncdataset.lat.data=curlat[:,thisy,thisx].reshape((ntimes,1,1))
            ncdataset.lon.data=curlon[:,thisy,thisx].reshape((ntimes,1,1))
            # set up this specific output filename
            filename=outputfile_base+"_sn{0:03}_we{1:03}".format(y,x)
            # now write the actual file. 
            write_file(filename,ncdataset,variable)
    if verbose:
        print("-------- Completed --------")

if __name__ == '__main__':
    try:
        # set up command line arguments and parse them
        parser= argparse.ArgumentParser(description="""Separates a 4D WRF grid (w/multiple files) into individual files
                                                       with timeseries, one for each x,y column. """)
        parser.add_argument('var',action='store',nargs="?",default="U",
                            help="Chose a variable variable to process [U,V]")
        parser.add_argument('nfiles',action='store',nargs="?",default="-1",
                            help="Number of files to process at a time, (<=0 means all files default=-1)")
        parser.add_argument('-v', '--version',action='version',version='wrf2anen.py 1.0')
        parser.add_argument ('--verbose', action='store_true',default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
        
        output_location="/glade/scratch/gutmann/swim/anen/extracting_"+str(res)+"km/processing/"+args.var+"/"
        nfiles=int(args.nfiles)
        if nfiles<=0:
            # process all WRF files at once (can take a long time)
            exit_code = main(args.var,output_location,filestart=0,fileend=None,verbose=args.verbose)
        else:
            # this hadn't really been debugged... it might work? 
            # process nfiles at a time still takes a while but memory requirements are lower and you can see output as you go
            # but you have to nccat the output files together at the end. 
            exit_code=1
            total_files=len(glob(wrf_location+wrfsearch))
            for i in range(0,total_files,nfiles):
                print(i)
                sys.stdout.flush()
                os.mkdir(output_location+str(i))
                if i+nfiles>total_files:
                    endfile=total_files
                else:
                    endfile=i+nfiles
                exit_code = main(args.var,output_location+str(i)+"/",filestart=i,fileend=endfile,verbose=args.verbose)
        
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
    