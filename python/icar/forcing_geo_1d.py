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
import datetime

import mygis

global verbose
verbose=False

geo_var_list=["lat","lon",'level']

def main (filename, outputfile):

    d = mygis.Dataset(filename)
    output = mygis.Dataset(outputfile,'w',format="NETCDF3_CLASSIC")

    for dim in d.dimensions:
        if dim=="time":
            output.createDimension(dim)
        else:
            output.createDimension(dim, d.dimensions[dim].size)


    for vname in d.variables:
        v=d.variables[vname]

        if not (v.name in geo_var_list):
            ov = output.createVariable(v.name, v.dtype, v.dimensions)
            ov[:] = v[:]
        else:
            if v.name == "lat":
                ov = output.createVariable("lat",v.dtype,(v.dimensions[0],))
                ov[:] = v[:,0]

            elif v.name == "lon":
                ov = output.createVariable("lon",v.dtype,(v.dimensions[1],))
                print("creating lon")
                print(v.dtype)
                print(v.dimensions[1])
                print(v[0,:])
                print(ov.shape)
                print(v.shape)
                print(ov)
                ov[:] = v[0,:]
            elif v.name == "level":
                ov = output.createVariable("z",v.dtype,v.dimensions)
                ov[:] = v[:]

        for a in v.ncattrs():
            try:
                ov.setncattr(a,v.getncattr(a))
            except Exception as e:
                print(e, a)

    for a in d.ncattrs():
        output.setncattr(a,d.getncattr(a))

    new_history = ', "'+" ".join(sys.argv)+'" at '+str(datetime.datetime.now())
    try:
        output.setncattr("history",output.getncattr("history")+new_history)
    except:
        output.setncattr("history",new_history)

    output.close()
    d.close



if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename',  action='store', default="some_file_name")
        parser.add_argument('outputfile',action='store', default="some_file_name")
        parser.add_argument('-v', '--version',action='version',
                version='Template Parser 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose=args.verbose

        verbose = args.verbose

        exit_code = main(args.filename, args.outputfile)
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
