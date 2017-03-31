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

global verbose
verbose=False

def main (filename):

    # TODO: Do something more interesting here...
    print('You want me to work on the file: '+filename)

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('filename',nargs="?", action='store', default="some_file_name")
        parser.add_argument('-f2', dest="file_two", action='store', default="a_Second_filename")
        parser.add_argument('-v', '--version',action='version',
                version='Template Parser 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose=args.verbose

        verbose = args.verbose

        exit_code = main(args.filename)
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
