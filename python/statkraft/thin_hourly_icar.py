#!/usr/bin/env python

"""
SYNOPSIS

    template_argparse.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION

    TODO This describes how to use this script.
    This docstring will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION


"""
from __future__ import absolute_import, print_function, division

import sys
import os
import traceback
import argparse
import glob

import numpy as np

import mygis
from bunch import Bunch

global verbose
verbose=False

icar_name_search = "icar_*2d.nc"

vars_to_copy = ["time","rain","snow","ta2m", "hus2m","u10m","rsds","rlds"]

# "lat","lon","time", will be copied explicitly
vars_to_write = ["rain","snow","ta2m", "hus2m","u10m","rsds","rlds"]
varname_translate=Bunch(u10m = "wind")

global_atts = Bunch(Conventions = "CF-1.6",
                    title = "Intermediate Complexity Atmospheric Research Model output",
                    institution = "National Center for Atmospheric Research",
                    source = "Intermediate Complexity Atmospheric Model version:0.9.1",
                    references = "Gutmann et al. 2016: The Intermediate Complexity Atmospheric Model. JHM doi:10.1175/JHM-D-15-0155.1",
                    contact = "Ethan Gutmann : gutmann@ucar.edu",
                    git = "master-30-g49f3d90")

def find_icar_files(dirname):
    files=glob.glob(dirname+"/statkraft/output/"+icar_name_search)
    files.sort()
    start_files=glob.glob(dirname+"/statkraft/output/start_jan_2020/"+icar_name_search)
    start_files.sort()
    for i in range(len(start_files)-1):
        curf=start_files[i].split("/")[-1]
        exists_in_files=False
        for j in range(len(files)):
            testf=files[j].split("/")[-1]
            if testf==curf:
                exists_in_files=True
                files[j] = start_files[i]
        if not exists_in_files:
            files.insert(i,start_files[i])

    if verbose:
        print("  Number of files:"+str(len(files)))
        print("  First file = "+files[0])
        print("  Last file = "+files[-1])
        print("  ")
    
    return files
        
def init_data(filename):
    lat = mygis.read_nc(filename,"lat")
    lon = mygis.read_nc(filename,"lon")
    
    last_rain = np.zeros(lat.data.shape)
    last_snow = np.zeros(lat.data.shape)
    
    return Bunch(last_rain=last_rain, last_snow=last_snow, lat=lat, lon=lon)

def update_snow(data):
    next_last_snow = np.zeros(data["lat"].data.shape)
    
    for i in range(data["snow"].data.shape[0]):
        next_last_snow[:] = data["snow"].data[i]
        cur_snow = data["snow"].data[i] - data.last_snow
        
        # if we have "negative" snow, then we probably just stepped back down to 0 accumulations
        if cur_snow.mean()<(-10): # use -10 to avoid possibly tripping over -1e-10 or the like
            if verbose:print("Found negative snow"+str(i))
            if (data["snow"].data[i].mean() < 100):
                cur_snow = data["snow"].data[i]

        # sanity check
        cur_snow[cur_snow<0] = 0
        data["snow"].data[i] = cur_snow
        
        #  this will leave it so that last_snow is always upto date. 
        data.last_snow[:] = next_last_snow

def update_rain(data, hourly_rain):
    next_last_rain = np.zeros(data["lat"].data.shape)
    
    for i in range(data["rain"].data.shape[0]):
        next_last_rain[:] = data["rain"].data[i]
        cur_rain = data["rain"].data[i] - data.last_rain
        
        # if we have "negative" rain, then we probably just stepped back down to 0 accumulations
        if cur_rain.mean()<(-10):
            if verbose:print("Found negative rain"+str(i))
            if (data["rain"].data[i].mean() < 100):
                cur_rain = data["rain"].data[i]
        elif hourly_rain[i].mean()>0.0001:
            cur_rain = hourly_rain[i]

        # sanity check
        cur_rain[cur_rain<0] = 0
        data["rain"].data[i] = cur_rain    
        
        #  this will leave it so that last_rain is always up to date. 
        data.last_rain[:] = next_last_rain


def update_data(filename, data):
    for v in vars_to_copy:
        if verbose:print("  "+v)
        data[v] = mygis.read_nc(filename,v)
        
        if v=="snow":
            update_snow(data)
        
        if v=="rain":
            update_rain(data, mygis.read_nc(filename,"rain_rate").data)
            
        if v=="u10m":
            data[v].data=np.sqrt(data[v].data**2 + mygis.read_nc(filename,"v10m").data**2)
            data[v].atts["standard_name"] = "10m_wind_speed"
            data[v].atts["long_name"] = "Wind speed at 10m"
        
def add_var(NCfile, varname, data, dtype, dims):
    NCfile.createVariable(varname,dtype, dims)
    NCfile.variables[varname][:] = data.data[:]
    for k in data.atts.keys():
        NCfile.variables[varname].__setattr__(k,data.atts[k])


def write_data(data, output_dir, filename):
    if verbose:print("Writing:"+filename)
    NCfile = mygis.Dataset(output_dir+"/"+filename, mode="w", format="NETCDF4")
    
    NCfile.createDimension("time", data["rain"].data.shape[0])
    NCfile.createDimension("lat", data["rain"].data.shape[1])
    NCfile.createDimension("lon", data["rain"].data.shape[2])
    
    varname="time"
    add_var(NCfile, varname, data[varname], "d", ("time",))

    varname="lat"
    add_var(NCfile, varname, data[varname], "f", ("lat","lon"))
    varname="lon"
    add_var(NCfile, varname, data[varname], "f", ("lat","lon"))

    dims=("time","lat","lon")
    dtype="f"
    for v in vars_to_write:
        try:
            outputvar = varname_translate[v]
        except:
            outputvar = v
        if verbose:print("  "+v+" = "+outputvar)
        add_var(NCfile, outputvar, data[v], dtype, dims)
    
    # for g in global_atts:
    
    NCfile.close()


def main (directory, output_dir):

    if verbose: print("Working on "+directory)
    
    files=find_icar_files(directory)
    
    data = init_data(files[0])
    
    for f in files:#[:2]:
        if verbose:print("Working on :"+f)
        update_data(f, data)
        write_data(data, output_dir, f.split("/")[-1])
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('directory',nargs="?", action='store', default="")
        parser.add_argument('-output', dest="output_dir", action='store', default=None)
        parser.add_argument('-v', '--version',action='version',
                version='Template Parser 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose = args.verbose
        
        if args.output_dir==None:
            output_dir = args.directory.split("/")[1]
        else:
            output_dir = args.output_dir
            
        if verbose:print("Outputing to:"+output_dir)
        try:
            os.makedirs(output_dir)
        except:
            pass

        exit_code = main(args.directory, output_dir)
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
