#!/usr/bin/env python
"""
SYNOPSIS

    mkscratch.py <directory_name>

DESCRIPTION

    Create and link to a directory in the yellowstone scratch filesystem
    If it doesn't exist already, creates a directory in scratch that matches
    the user's current directory in project space plus the new directory.
    Then create a symlink to the directory_name specified on the commandline.
    in the current directory.

EXAMPLES

    mkscratch.py temp_output

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


def main(newdir, verbose, quiet):
    """docstring for main"""

    if verbose:print("---------------------------")
    # get the current working directory, resolving symlinks as necessary
    working_dir=os.path.realpath(os.getcwd())

    scratch_dir = working_dir.replace("/glade/p/work/","/glade2/scratch2/")
    scratch_dir = scratch_dir.replace("/glade2/work/","/glade2/scratch2/") # updated work directory
    scratch_dir+="/"+newdir

    if verbose:
        print("Currently in : "+working_dir)
        print("Creating link to : "+scratch_dir)

    # Error checking, does the directory exist already? Worse, is it a file?
    if os.path.isdir(scratch_dir):
        if not quiet:
            print("  WARNING: Scratch directory already exists!")
            print("    Creating link, but be careful not to over write existing files. ")
    if os.path.isfile(scratch_dir):
        if not quiet:
            print("ERROR: Scratch name is an existing file!")
            print(scratch_dir)
        raise SystemExit("Invalid Directory Name")

    fullpath=""
    # work through the path creating directories as necessary
    for curdir in scratch_dir.split("/"):
        # add the next level in the directory hierarchy to check
        fullpath+="/"+curdir

        # Error checking, is it an existing file?
        if os.path.isfile(fullpath):
            if not quiet:
                print("ERROR: Scratch name is an existing file!")
                print(fullpath)
            raise SystemExit("Invalid Directory Name")

        # does this directory exist already?  If not, create it.
        if not os.path.isdir(fullpath):
            if verbose:print("Creating Scratch Directory:"+curdir)
            os.mkdir(fullpath)

    if verbose:
        print("Setting up link for:"+newdir)
        print("---------------------------")

    os.symlink(scratch_dir,newdir)



if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Create and link to a scratch dir',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename', action='store')
        parser.add_argument('-v', '--version',action='version',
                version='mkscratch 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        parser.add_argument ('-q','--quiet', action='store_true',
                default=False, help='quiet output', dest='quiet')
        args = parser.parse_args()

        exit_code = main(args.filename, args.verbose, args.quiet)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print("---------------------------")
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        print("---------------------------")
        traceback.print_exc()
        os._exit(1)
