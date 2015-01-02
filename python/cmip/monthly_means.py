#!/usr/bin/env python
from __future__ import print_function
import glob
import sys
import mygis


marchstart=(31+28)*24
marchend=(31+28+31)*24

startdate,enddate=marchstart,marchend
outputfilename="mean_{}.nc"

def load_temperature(filename):
    """docstring for load_temperature"""
    pot_t=mygis.read_nc(filename,"th",returnNCvar=True)
    pii=mygis.read_nc(filename,"pii",returnNCvar=True)
    tskin=mygis.read_nc(filename,"tskin").data
    
    tskin+=pot_t.data[:,0,:]*pii.data[:,0,:]
    pot_t.ncfile.close()
    pii.ncfile.close()
    return (tskin/2)-273.15

def load_precip(filename):
    """docstring for load_precip"""
    cur_data=mygis.read_nc(filename,"rain").data
    step=int(filename.split("-")[-1])
    if step==0:
        return cur_data
    else:
        last_file="-".join([filename.split("-")[:-1], "{:05}".format(step-1)])
        last_data=mygis.read_nc(last_file,"rain").data
        if cur_data.sum()<last_data.sum():
            return cur_data
        else:
            return cur_data-last_data


def load_var(filename,varname):
    """docstring for load_var"""
    if varname=="t":
        return load_temperature(filename)
    elif varname=="rain":
        return load_precip(filename)
    else:
        return mygis.read_nc(filename,varname).data

def main():
    """docstring for main"""
    files=glob.glob("output/icar_*-0*")
    files.sort()
    means=dict()
    variables=["t","rain"]
    n=0
    for f in files:
        step=int(f.split("-")[-1])
        if (step>=startdate) and (step<=enddate):
            n+=1
            print(f,end="\r")
            sys.stdout.flush()
            for curvar in variables:
                means[curvar]+=load_var(f,curvar)
    print("")
    print("Number of data points: {}".format(n))
    for curvar in variables:
        means[curvar]/=n
        mygis.write(outputfilename.format(curvar),means[curvar])
    
        

if __name__ == '__main__':
    main()