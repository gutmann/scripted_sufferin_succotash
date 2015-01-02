#!/usr/bin/env python
# encoding: utf-8
"""
snow_raw2data.py

USAGE: 
snow_raw2data [-o outputfile] raw_datafile keyfile
    raw_datafile: the intput data file to be processed
        This file must be a whitespace delimited column format file with no missing cells
        (missing cells should be filled with e.g. -9999)
    keyfile: a file that contains a description of the raw_datafile as follows
        line1: a number specifying the date format for this file
            1 = Modified Julian Day
            2 = Julian Day
            3 = Excel Date time stamp PC format
            4 = Excel Date time stamp Mac format
            5 = 5 column Year, Month, Day, Hour, Minute
            6 = 3 column Year, Day-of-year, HHmm (e.g. 10:30AM would be 1030)
            if the number is negative it means the date format is stored in another file 
            (e.g. when the data come from a CR10x with a separate table)
        if line1 is negative: 
            line2: the name of the date file
        else: 
            line2 is ignored
        
        line3: the column in the date file for the specified date (column numbers are zero based)
            if the specified date is Y,M,D,hr,mn then column numbers must be listed in that order
        line4: columns of data to be read and output (column numbers are zero based)
        line5: an offset to the time in days

    -o outputfile: optionally an output filename may be specified.  If none is specified
        the default is raw_datafile.out

Created by Ethan Gutmann on 2011-04-20.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""

import sys
import getopt
import numpy as np
import load_data
import date_fun
from bunch import Bunch

#this should be in the DOCSTRING above instead of here... eventually 
# transitioning to a new template file will fix that...
help_message = '''
snow_raw2data [-o outputfile] raw_datafile keyfile
    raw_datafile: the intput data file to be processed
        This file must be a whitespace delimited column format file with no missing cells
        (missing cells should be filled with e.g. -9999)
    keyfile: a file that contains a description of the raw_datafile as follows
        line1: a number specifying the date format for this file
            1 = Modified Julian Day
            2 = Julian Day
            3 = Excel Date time stamp PC format
            4 = Excel Date time stamp Mac format
            5 = 5 column Year, Month, Day, Hour, Minute
            6 = 3 column Year, Day-of-year, HHmm (e.g. 10:30AM would be 1030)
            if the number is negative it means the date format is stored in another file 
            (e.g. when the data come from a CR10x with a separate table)
        if line1 is negative: 
            line2: the name of the date file
        else: 
            line2 is ignored
        
        line3: the column in the date file for the specified date (column numbers are zero based)
            if the specified date is Y,M,D,hr,mn then column numbers must be listed in that order
        line4: columns of data to be read and output (column numbers are zero based)
        line5: an offset to the time in days

    -o outputfile: optionally an output filename may be specified.  If none is specified
        the default is raw_datafile.out
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


# parse_keys(keyfile) takes a filename and reads the keys as defined in the help message
def parse_keys(keyfile):
    f=open(keyfile,'r')
    # read each line into it's own variable
    l1=f.next()
    l2=f.next()
    l3=f.next()
    l4=f.next()
    try: # because the time offset was added after this program had been used a bit we catch lots of exceptions
        l5=np.float(f.next().split()[0])
    except ValueError: # if this line does not contain a number np.float will fail
        l5=None
    except StopIteration: # if this line is past the end of the file f.next will fail
        l5=None
    except IndexError: # if this line is blank split()[0] will fail
        l5=None
    try: # also, the time limits were added after this program had been used a bit we catch lots of exceptions
        l6=np.array(f.next().split()).astype("f")
        limit1=date_fun.date2mjd(l6[0],l6[1],l6[2],0,0,0)
        limit2=date_fun.date2mjd(l6[3],l6[4],l6[5],0,0,0)
        l6=np.array([limit1,limit2])
    except ValueError: # if this line does not contain a number astype will fail?
        l6=[0]
    except StopIteration: # if this line is past the end of the file f.next will fail
        l6=[0]
    except IndexError: # if this line is blank split()[0] will fail?
        l6=[0]
    # close the input file
    f.close()

    # line 1 is the datetype, 1=MJD,2=JD,3=Excel_pc,4=Excel_mac,5=Y M D h m
    datetype=int(l1.split()[0])
    # negative datetype values mean that the dates come from a separate file
    # to be read from line 2
    if datetype<0 : 
        datefile=l2.split()[0]
    else: 
        datefile=''
    # line 3 stores the column numbers for the date column (one column for datetype=1-4, 5 for datetype=5,3 for dt=6)
    datecols=np.array(l3.split()).astype('i')
    datacols=np.array(l4.split()).astype('i')
    time_offset=l5
    timelimits=l6
    
    # return what will look like a C-struct, and remove the negative sign from datetype if necessary
    return Bunch(datetype=abs(datetype),datefile=datefile,datecols=datecols,datacols=datacols,
                 time_offset=time_offset,timelimits=timelimits)

# read the dates from a set of data using the datetype and columns specified by the keys
def readdates(data,keys):
    dates=data[:,keys.datecols]
    # we don't need to do anything if datetype=1 (we want to return Modified Julian Days)
    if keys.datetype==1 :
        outputdates=dates
    # if dates are Julian Days, subtract 2400000.5 to make them modified Julian days
    if keys.datetype==2 :
        outputdates=dates-2400000.5
    # if dates are excel PC format (days from 1899 Dec.31?), use date_fun to convert them modified Julian days
    elif keys.datetype==3 :
        outputdates=date_fun.excel2mjd(dates)
    # if dates are excel Mac format (days from 1904 Jan.1), use date_fun to convert them modified Julian days
    elif keys.datetype==4 :
        outputdates=date_fun.excel2mjd(dates)
    # if dates are Year, Month, Day, Hour, Minute use date_fun to conver them to modified Julian days
    elif keys.datetype==5 :
        yearoffset=0
        if np.max(dates[:,0])<1000:yearoffset=2000
        outputdates=date_fun.date2mjd(dates[:,0]+yearoffset,dates[:,1],dates[:,2],dates[:,3],dates[:,4])
    elif keys.datetype==6 :
        outputdates=date_fun.ydoyhm2mjd(dates[:,0],dates[:,1],dates[:,2])
    else : 
        print("ERROR: unknown date type")
        raise(ValueError)
    outputdates.shape=(len(outputdates),1)
    return outputdates

        

def snow_raw2data(inputfile,keyfile,outputfile):
    # read information about the input file from the keyfile
    keys=parse_keys(keyfile)
    # if dates are stored in a separate datefile, read that first
    if keys.datefile:
        # load dates
        dates=load_data.cols(keys.datefile,dtype='d')
        # then load data
        data=load_data.cols(inputfile)
        # read the dates from the datefile
        dates=readdates(dates,keys)+keys.time_offset
    else:
        # load data, dates are stored within data
        data=load_data.cols(inputfile,dtype='d')
        # read the dates from the datafile
        dates=readdates(data,keys)
    if len(keys.timelimits)>1:
        # find the dates in between timelimits[0] and timelimits[1]
        # index the where tuple along the 0th axis and at the first point in the array
        startpt=np.where(dates>=keys.timelimits[0])[0][0]
        stoppt=np.where(dates>=keys.timelimits[1])[0][0]
        data=data[startpt:stoppt,:]
        dates=dates[startpt:stoppt]
    
    data=data[:,keys.datacols]
    outputdata=np.concatenate((dates,data),axis=1)
    
    np.savetxt(outputfile,outputdata,fmt='%10.5f')

    

# main program, just parses command line arguments
def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
    
        # option processing
        output=None
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
        if len(args)==0:
            raise Usage(help_message)
    
    except Usage, err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg))
        print("\t for help use --help")
        return 2
    
    # if re.match(args[0],".*\.dat"):
    #     pass #run csi2space on data file?
    if output==None :
        output=args[0]+'.out'
    # call processing program
    snow_raw2data(args[0],args[1],output)
    # try:
    #   snow_raw2data(args[0],args[1],output)
    # except ValueError:
    #   print("quitting...")
        




if __name__ == "__main__":
    sys.exit(main())
