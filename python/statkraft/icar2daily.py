#!/usr/bin/env python

"""
SYNOPSIS

    icar2daily.py [-h] [--verbose] [-v, --version] <filename> [-o outputfile]

DESCRIPTION

    Reads one or more ICAR output files and creates daily precip, tmin, and tmax files
    
EXAMPLES

    icar2daily.py output/icar_2021_\*  -o rcp45_2021

EXIT STATUS

    0=success
    1=error

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    1.0
    
"""

import sys
import os
import traceback
import argparse
import glob

import numpy as np
import mygis
from bunch import Bunch

import date_fun

global verbose
verbose=True

missing = -999.99
mask_file = "geo_em.d02.nc"

prec_atts = Bunch(temporal_period="00-24Z", missing_value=missing, _FillValue=missing, description="Daily precipitation amount", units="mm")
tmin_atts = Bunch(temporal_period="00-24Z", missing_value=missing, _FillValue=missing, description="Daily minimum temperature",  units="degrees C")
tmax_atts = Bunch(temporal_period="00-24Z", missing_value=missing, _FillValue=missing, description="Daily maximum temperature",  units="degrees C")
tave_atts = Bunch(temporal_period="00-24Z", missing_value=missing, _FillValue=missing, description="Daily average temperature",  units="degrees C")

lat_atts= Bunch(axis="Y",datum="WGS84", units="degrees_north", long_name="Latitude", standard_name="latitude_north")
latvar  = Bunch(data=None, attributes=lat_atts, dims=("lat","lon"), dtype='f', name="lat")

lon_atts= Bunch(axis="X",datum="WGS84", units="degrees_east", long_name="Longitude", standard_name="longitude_east")
lonvar  = Bunch(data=None, attributes=lon_atts, dims=("lat","lon"), dtype='f', name="lon")

time_atts= Bunch(axis="T", units="days since 1800-01-01 00:00:00", calendar="noleap", long_name="time", standard_name="time")
era_time_atts= Bunch(axis="T", units="days since 1858-11-17 00:00:00", calendar="gregorian", long_name="modified Julian Day", standard_name="time",UTCoffset=0)
timevar  = Bunch(data=None, attributes=time_atts, dims=("time",), dtype='f', name="time")

space_time_vars = [latvar, lonvar, timevar]

global_attributes=Bunch(institute = "National Center for Atmospheric Research", 
                        contact   = "Ethan Gutmann: gutmann@ucar.edu", 
                        source    = "The Intermediate Complexity Atmospheric Research (ICAR) model",
                        references= "Gutmann et al. 2016: The Intermediate Complexity Atmospheric Research Model. J.Hydrometeor. doi:10.1175/JHM-D-15-0155.1, in press.",
                        Conventions = "CF-1.6", 
                        scenario  = "",
                        forcing   = "")

class ICAR_Reader(object):
    files=[]
    curstep=0
    curfile=0
    rr_var="rain_rate"
    ac_var="rain"
    t_var="ta2m"
    
    def __init__(self, filesearch, steps_per_day=24):
        """docstring for init"""
        self.files=glob.glob(filesearch)
        self.files.sort()
        self.steps_per_day=steps_per_day
        
        if len(self.files)==0:
            raise ValueError("\nERROR: no files match search string:"+filesearch+"\n")
        if verbose: print("Reading:{}".format(self.files[self.curfile]))
        self.tdata  = mygis.read_nc(self.files[self.curfile], self.t_var).data
        # self.rrdata = mygis.read_nc(self.files[self.curfile],self.rr_var).data
        # self.acdata = mygis.read_nc(self.files[self.curfile],self.ac_var).data
        
        # self.mask = mygis.read_nc(mask_file,"LANDMASK").data[0]==0
        
        
    def update_data(self):
        """docstring for update_data"""
        if verbose: print("Finished:{}".format(self.files[self.curfile]))
        self.curfile+=1
        self.curstep=0
        
        if self.curfile==len(self.files):
            raise StopIteration
            
        if verbose: print("Reading:{}".format(self.files[self.curfile]))
        next_tdata  = mygis.read_nc(self.files[self.curfile], self.t_var).data
        # next_rrdata = mygis.read_nc(self.files[self.curfile],self.rr_var).data
        # next_acdata = mygis.read_nc(self.files[self.curfile],self.ac_var).data
        # check that the next file has at least one day's worth of data in it
        # if next_rrdata.shape[0]<24:
        #     raise StopIteration
            
        # fill in any missing dates in the rain_rate data (happens with restarts)
        # for i in range(next_rrdata.shape[0]):
        #     if next_rrdata[i].max()==0:
        #         if i==0:
        #             next_rrdata[i]=next_acdata[i]-self.acdata[-1]
        #         else:
        #             next_rrdata[i]=next_acdata[i]-next_acdata[i-1]
        #         if verbose:
        #             if next_rrdata[i].max()>0:
        #                 print("Found a restart step at step {} in file:{}".format(i,self.files[self.curfile]))            
                        
        self.tdata  = np.ma.array(next_tdata, mask=next_tdata>10000)
        # self.acdata = np.ma.array(next_acdata, mask=next_acdata>10000)
        # self.rrdata = np.ma.array(next_rrdata, mask=next_rrdata>10000)
        
        
    def next(self):
        """docstring for next"""
        if self.curstep==self.tdata.shape[0]:
            self.update_data()
        
        if self.curstep+self.steps_per_day > self.tdata.shape[0]:
            print("Ran out of files to process at step {} in file:{}".format(self.curstep,self.files[self.curfile]))
            print("If isn't the last file ({}), may need to edit to permit fractional days / file".format(self.files[-1]))
            raise StopIteration

        rain=[]
        tmin=[]
        tmax=[]
        # rain = self.rrdata[self.curstep : self.curstep+self.steps_per_day].sum(axis=0)
        # tmin = self.tdata[ self.curstep : self.curstep+self.steps_per_day].min(axis=0)  - 273.15
        # tmax = self.tdata[ self.curstep : self.curstep+self.steps_per_day].max(axis=0)  - 273.15
        tave = self.tdata[ self.curstep : self.curstep+self.steps_per_day].mean(axis=0) #- 273.15
        
        self.curstep+=self.steps_per_day
        
        # tmin[self.mask]=missing
        # tmax[self.mask]=missing
        # tave[self.mask]=missing
        
        # rain[rain>3000] = missing
        # tmin[tmin>100]  = missing
        # tmax[tmax>100]  = missing
        tave[tave>100]  = missing
        
        return (rain, tmin, tmax, tave)
    
    def __iter__(self):
        """docstring for __iter__"""
        return self
        
    # python3x compatibility
    __next__=next
    

def write_file(fname,data, varname, varatts):
    """docstring for write_file"""
    outputdata=np.zeros((len(data),data[0].shape[0],data[0].shape[1]))
    for i in range(len(data)):
        outputdata[i]=data[i]
    mygis.write(fname,outputdata, dtype='f', dims=("time","lat","lon"), varname=varname, attributes=varatts,
                extravars=space_time_vars, global_attributes=global_attributes)
    

def main (filename, outputfile):
    
    icar_data=ICAR_Reader(filename)
    if verbose:print(icar_data.files[0],icar_data.files[-1])
    
    raindata=[]
    tmindata=[]
    tmaxdata=[]
    tavedata=[]
    if verbose:print("Looping through data")
    for data in icar_data:
        raindata.append(data[0])
        tmindata.append(data[1])
        tmaxdata.append(data[2])
        tavedata.append(data[3])
    
    latvar.data=mygis.read_nc(icar_data.files[0],"lat").data
    lonvar.data=mygis.read_nc(icar_data.files[0],"lon").data
    
    dates_are_mjd = (icar_data.files[0].split("/")[1] == "erai")
    if verbose:print("using modified julian day dates="+str(dates_are_mjd))
    year = int(icar_data.files[0].split("/")[-1].split("_")[1])
    if dates_are_mjd:
        timevar.attributes = era_time_atts
        start_date = date_fun.date2mjd(year,01,01,00,00)
    else:
        start_date = (year-1800) * 365
    timevar.data = np.arange(start_date,start_date+len(raindata))
    if verbose:print(year, timevar.data[0], len(raindata))
    
    if verbose:print("Writing data")
    # write_file(outputfile+"_rain.nc",raindata, varname="precipitation_amount",      varatts=prec_atts)
    # write_file(outputfile+"_tmin.nc",tmindata, varname="daily_minimum_temperature", varatts=tmin_atts)
    # write_file(outputfile+"_tmax.nc",tmaxdata, varname="daily_maximum_temperature", varatts=tmax_atts)
    write_file(outputfile+"_tave",tavedata, varname="daily_average_temperature", varatts=tave_atts)
    
    
    
if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Generate daily rain, tmin, tmax, and tave files from hourly ICAR output. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename',action='store', help="glob file search for input files")
        parser.add_argument('-o', dest="outputfile", action='store', default="daily_", help="outputfile prefix")
        parser.add_argument('-m', dest="model", action='store', default="", help="Model used to drive ICAR, e.g. icar_run_cesm_rcp45")
        parser.add_argument('-v', '--version',action='version',
                version='icar2daily 1.1')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
        verbose=args.verbose
        
        print(args.model)
        try:
            global_attributes.scenario = args.model.split("_")[3]
        except:
            #era data
            global_attributes.scenario="Historical"
        
        if global_attributes.scenario[0]!="r":global_attributes.scenario="Historical"
        global_attributes.forcing  = args.model.split("_")[2]
        
        exit_code = main(args.filename,args.outputfile)
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
