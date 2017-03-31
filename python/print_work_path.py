#!/usr/bin/env python
"""
SYNOPSIS

    print_work_path.py

DESCRIPTION

    Print the path to a directory in the yellowstone p/work filesystem corresponding to the current scratch path
    If it doesn't exist already, creates a directory in work that matches
    the user's current directory in scratch space. 

EXAMPLES

    cd `print_work_path.py`

EXIT STATUS

    None

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


def main(verbose):
    """docstring for main"""

    if verbose:print("---------------------------")
    # get the current working directory, resolving symlinks as necessary
    scratch_dir=os.path.realpath(os.getcwd())

    work_dir = scratch_dir.replace("scratch","p/work")

    if verbose:
        print("Currently in : "+scratch_dir)
        print("Creating printing path to : "+work_dir)

    # Error checking, does the directory exist already? Worse, is it a file?
    if verbose:
        if os.path.isdir(work_dir):
            print("  Work directory already exists!")
    if os.path.isfile(work_dir):
        print("ERROR: Work name is an existing file!")
        print(work_dir)
        raise SystemExit("Invalid Directory")
        
    fullpath=""
    # work through the path creating directories as necessary
    for curdir in work_dir.split("/"):
        # add the next level in the directory hierarchy to check
        fullpath+="/"+curdir
        
        # Error checking, is it an existing file?
        if os.path.isfile(fullpath):
            print("ERROR: Work name is an existing file!")
            print(fullpath)
            raise SystemExit("Invalid Directory Name")
        
        # does this directory exist already?  If not, create it. 
        if not os.path.isdir(fullpath):
            if verbose:print("Creating Work Directory:"+curdir)
            os.mkdir(fullpath)

    print(work_dir)


if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Create and print path to a work dir',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-v', '--version',action='version',
                version='print_work_path 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        exit_code = main(args.verbose)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print("---------------------------")
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        print("---------------------------")
        traceback.print_exc()
        os._exit(1)
