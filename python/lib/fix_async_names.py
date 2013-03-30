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
import glob
import os

import numpy as np

filesearch="*.nc"
addstring="4km"
position=1
striplocations=2
def main (directory,prefix=""):
    curdir=os.getcwd()
    os.chdir(directory)
    files=glob.glob(filesearch)
    for f in files:
        fileparts=f.split("_")[striplocations:]
        newfile=prefix+"_".join(np.hstack([fileparts[:position],[addstring],fileparts[position:]]))
        print(f,newfile)
        os.rename(f,newfile)
        
    os.chdir(curdir)
        

    
    
if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Batch process filenames from asyncronous regression. ')
        parser.add_argument('directory',action='store')
        parser.add_argument('-v', '--version',action='version',
                version='fix_async_names 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        exit_code = main(args.directory)
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
