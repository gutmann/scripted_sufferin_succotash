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

import sys
import os
import traceback
import argparse

import mygis
import matplotlib.pyplot as plt

def main (filename, varname, d0, d2, outputfile, vmin,vmax):

    data=mygis.read_nc(filename,varname).data
    if d2!=None:
        data=data[:,:,int(d2)]
    if d0!=None:
        data=data[int(d0)]

    plt.imshow(data,vmin=vmin,vmax=vmax)
    plt.colorbar()
        
    plt.savefig(outputfile)

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename', action='store', default="some_file_name")
        parser.add_argument('varname',  action='store', default="data")
        parser.add_argument('-d0',nargs="?", action='store', default=None)
        parser.add_argument('-d2',nargs="?", action='store', default=None)
        parser.add_argument('-vmin',nargs="?", action='store', default=None,type=int)
        parser.add_argument('-vmax',nargs="?", action='store', default=None,type=int)
        parser.add_argument('-o', nargs="?", action='store', default="output.png")
        parser.add_argument('-v', '--version',action='version',
                version='Template Parser 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        exit_code = main(args.filename, args.varname,args.d0, args.d2, args.o, args.vmin, args.vmax)
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
