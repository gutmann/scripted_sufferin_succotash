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

import xarray as xr
import numpy as np

from scipy.interpolate import interp2d, interp1d

global verbose
verbose=False

vars_to_correct = ["qv"]#, "th", "p", "u", "v"]
vars_to_correct = ["qv","u","v","theta","p","swdown","lwdown","tskin"]
z_var_name = "z"

def z_interp(inputdata, outputdata, lut):

    if (len(outputdata.shape)==4):
        for i in range(outputdata.shape[2]):
            for j in range(outputdata.shape[3]):
                outputdata[:,:,i,j] = ( (inputdata[:,lut[0,:,i,j].astype('i'), i,j]      * lut[np.newaxis, 1,:,i,j])
                                      + (inputdata[:,lut[0,:,i,j].astype('i')+1, i,j]) * (1 - lut[np.newaxis, 1,:,i,j]))
    elif (len(outputdata.shape)==3):
        for i in range(outputdata.shape[1]):
            for j in range(outputdata.shape[2]):
                outputdata[:,i,j] = ( (inputdata[lut[0,:,i,j].astype('i'), i,j]    * lut[1,:,i,j])
                                    + (inputdata[lut[0,:,i,j].astype('i')+1, i,j]) * (1 - lut[1,:,i,j]))

def compute_column_vlut(zin, zout):
    vlut=np.zeros((2,zout.shape[0]))
    last = 0

    for i in range(zout.shape[0]):
        while ((zin[last] < zout[i]) and (last < (zin.size)-1)):
            last+=1

        if zin[last] >= zout[i]:
            if (last==0):
                vlut[0,i] = 0
                vlut[1,i] = 1
            else:
                vlut[0,i] = last-1

                dz = zin[last]-zin[last-1]
                if dz == 0:
                    # should never get here, but good to check
                    vlut[0,i]=1
                else:
                    distance = zin[last]-zout[i]
                    vlut[1,i] = distance / dz
        else:
            if last != (zin.size-1):
                raise ValueError
            else:
                vlut[0,i] = -2
                vlut[1,i] = 0
    return vlut



def compute_vinterp(zin, zout):
    vlut = np.zeros((2,zout.shape[0],zout.shape[1],zout.shape[2]))
    print("Computing VLUT")
    zin.shape
    for i in range(zin.shape[1]):
        print(i, zin.shape[1])
        for j in range(zin.shape[2]):
            tmp = compute_column_vlut(zin[:,i,j], zout[:,i,j])
            vlut[:,:,i,j] = tmp

    return vlut

def z_interp_all_vars(base_data, output_data, zin, zout, varlist):

    vlut = compute_vinterp(zin, zout)

    for v in varlist:
        print("interpolating "+v)
        z_interp(base_data[v], output_data[v], vlut)



def main (era_filesearch, cesm_base_filesearch, bias_output):

    print("opening data")
    era_data         = xr.open_mfdataset(era_filesearch,         concat_dim='time')
    base_cesm_data   = xr.open_mfdataset(cesm_base_filesearch,   concat_dim='time')

    print("loading data")
    era_data.load()
    base_cesm_data.load()

    print("compute means")
    emean = era_data.std(dim="time")
    cmean = base_cesm_data.std(dim="time")

    print("creating data")
    interpolated_era = xr.zeros_like(cmean)
    print("loading data")
    interpolated_era.load()

    z_interp_all_vars(emean, interpolated_era, era_data["z"].mean(dim="time"), base_cesm_data["z"].mean(dim="time"), vars_to_correct)
    interpolated_era.to_netcdf("era_interpolated_std.nc")

    print("Computing Bias")
    bias = interpolated_era - cmean

    print("writing")
    bias.to_netcdf(bias_output)


if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # parser.add_argument('filename',nargs="?", action='store', default="some_file_name")
        parser.add_argument('-f1', dest="f1", action='store', default="")
        parser.add_argument('-f2', dest="f2", action='store', default="")
        parser.add_argument('-o',  dest="output", action='store', default="")
        parser.add_argument('-v', '--version',action='version',
                version='Template Parser 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose=args.verbose

        verbose = args.verbose

        exit_code = main(args.f1, args.f2, args.output)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
