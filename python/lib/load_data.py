'''
Created on Mar 18, 2011

@author: gutmann

    This module is for programs used to load data from various file formats

'''

import csv
import numpy as np
from datetime import datetime
import re
import subprocess
import shutil
import os
# import excel2csv
import glob

def cols_date(filename,dtype='f',year=-1,month=-1,day=-1,hour=-1,minute=-1,second=-1):
    '''
    cols: reads a text file with an arbitrary number of header lines
    starts at first line that has nothing but numbers on it. 
    returns a numpy array of data using numpy.loadtxt
    
    defaults to float data type, can be made to use other datatypes with dtype='d'
    
    EXAMPLE:
    import load_data
    data=load_data.cols('filename')
    '''
    from julday import mjul_day

    data=cols(filename,dtype=dtype)
    years=data[:,year]
    months=data[:,month]
    days=data[:,day]
    if hour>-1:
        hours = data[:,hour] 
    else: hours=np.zeros(len(days))
    if minute>-1:
        minutes = data[:,minute]
    else: minutes=np.zeros(len(days))
    if second>-1:
        seconds = data[:,second]
    else: seconds=np.zeros(len(days))
    times=[mjul_day(years[i], months[i], days[i], hours[i], minutes[i], seconds[i]) for i in range(len(years))]
    #    times=[datetime(result[:,year],result[:,month],result[:,day],result[:,hour],result[:,minute])]
    times=np.array(times)
    data=data[:,max((year,month,day,hour,minute))+1:]
    return (times, data)
#end


from bunch import Bunch
def cols(filename,dtype='d',delimiter=None,readheader=False):
    '''
    cols: reads a text file with an arbitrary number of header lines
    starts at first line that has nothing but numbers on it. 
    returns a numpy array of data using numpy.loadtxt
    
    defaults to float data type, can be made to use other datatypes with dtype='d'
    
    EXAMPLE:
    import load_data
    data=load_data.cols('filename')
    '''
    from numpy import loadtxt
    from is_number import is_number
    
    headerdata=''
    f=open(filename, 'r')
    
    inheader=True
    headerlength=-1
    while inheader:
        headerlength+=1
        line=f.readline()
        curdata=line.split()
        inheader=False
        for test in curdata:
            inheader= (inheader or (not is_number(test)))
        #endfor
        if inheader & readheader:
            headerdata+=line
    #endwhile
    
    f.close()
    
    if readheader:
        return Bunch(data=loadtxt(filename,skiprows=headerlength,dtype=dtype,delimiter=delimiter),header=headerdata)
    else:
        return loadtxt(filename,skiprows=headerlength,dtype=dtype,delimiter=delimiter)
#end


def readfirst(data):
    line=data.next()
    for i in range(len(line)):
        line[i]=re.sub('(?P<num>[0-9])-','\g<num> ',line[i])
        line[i]=line[i].replace('/',' ').replace(':',' ').replace('+',' ').replace('kg',' ').replace('\'',' ').replace('"',' ')
    time=np.array(line[0].split(),'i')
    if time[0]<1800: time[0]+=2000
    time=datetime(time[0],time[1],time[2],time[3],time[4])
    outputtime=[datetime(2000,1,1) for i in range(100)]
    outputtime[0]=time

    outputdata=np.zeros((100,len(line)-1),'d')
    outputdata[0,:]=np.array(line[1:],'d')
    return (outputtime,outputdata)
#end readfirst


def datalogger(filename, dtype='d', fill=-9999):
    '''
    Reads in a Campbell scientific datalogger file and returns a list
        (list of datetimes, array of data)
    by default returns doubles could add dtype='d', ...
    fill=value to fill empty cells ('') with default=-9999
        
    EXAMPLE: 
    import load_data
    (dates,data)=load_data.datalogger('CR1000_Table1.dat',fill=-9999.99)
    '''

    data=csv.reader(open(filename,'rU'))#,dialect=csv.excel)
    for i in range(4):line=data.next() # skip the four standard header lines

    (outputtime,outputdata)=readfirst(data)
    
    i=1l
    n=len(outputtime)
    ncols=len(outputdata[0,:])
    for line in data:
        if line and line[0]:
            for i in range(len(line)):
                line[i]=re.sub('(?P<num>[0-9])-','\g<num> ',line[i])
                line[i]=line[i].replace('/',' ').replace(':',' ').replace('+',' ').replace('kg',' ').replace('\'',' ').replace('"',' ')
            time=np.array(line[0].split(),'i')
            if time[0]<1800: time[0]+=2000
            outputtime[i]=datetime(time[0],time[1],time[2],time[3],time[4])
    
    
            try:
                outputdata[i,:]=np.array(line[1:ncols+1],dtype)
            except ValueError:
                outputdata[i,:]= -9999.9
    
            i+=1
            if (i>=n):
                newoutputdata=np.zeros((n*2,len(line)-1),dtype)
                newoutputdata[0:i]=outputdata
                newoutputdata[0:i]=outputdata
                outputdata=newoutputdata
    
                newoutputtime=[datetime(2000,1,1) for j in range(n*2)]
                newoutputtime[0:i]=outputtime
                outputtime=newoutputtime
                n*=2
            #endif i>=n
        #endif line[0]
    #endfor
    i-=1
    return (outputtime[0:i],outputdata[0:i,:])
#end load_datalogger

# def extractcsv(filename):
#     # initcsvs=glob.glob('*.csv')
#     excel2csv.excel2csv(filename,onesheet=1,outputfile='temporary.csv')
#     # newcsvs=glob.glob('*.csv')
#     # for curfile in newcsvs:
#     #     if curfile not in initcsvs:
#     #         return curfile
#     # return ''
#     return "temporary.csv"
    
    
def decagon(filename,fromexcel2csv=True):
    removefile=False
    if fromexcel2csv==False: 
        print("I'm sorry I don't know how to read in AM/PM date strings yet please fix me")
        return
    if re.match('.*\.xls$',filename):
        removefile=True
        filename=extractcsv(filename)
        if filename=='':
            print("Could not find output csv file!")
            return
    if re.match('.*\s',filename):
        newfile=filename.replace(' ','_')
        shutil.copyfile(filename,newfile)
        if removefile:
            os.remove(filename)
        filename=newfile
        removefile=True
    # removes commas from within the file and replace them with spaces works even if there are no commas to remove
    p=subprocess.Popen("comma2space "+filename+" "+filename+".txt",shell=True)
    p.wait()
    if removefile:
        os.remove(filename)
    filename=filename+".txt"
    removefile=True
    
    d=cols(filename)
    
    if removefile:
        os.remove(filename)
    return d

