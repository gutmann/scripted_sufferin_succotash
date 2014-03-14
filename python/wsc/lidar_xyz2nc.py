#!/usr/bin/env python
import glob,os
import numpy as np
import swim_io as io
from bunch import Bunch

dx=10 #m

def load_data(filename):
    try:
        data=np.loadtxt(filename,delimiter=",")
    except:
        data=np.loadtxt(filename)
    x=data[:,0]
    y=data[:,1]
    z=data[:,2]
    
    xmin=x.min()
    ymin=y.min()
    xsize=(x.max()-xmin) / dx+1
    ysize=(y.max()-ymin) / dx+1

    acc_snowdepth=np.zeros((ysize,xsize))
    n_observations=np.zeros((ysize,xsize))
    
    xmap=np.arange(xmin+dx/2,x.max()+dx/2,dx)
    ymap=np.arange(ymin+dx/2,y.max()+dx/2,dx)
    
    for i in range(len(x)):
        xpos=(x[i]-xmin)/dx
        ypos=(y[i]-ymin)/dx
        acc_snowdepth[ypos,xpos]+=z[i]
        n_observations[ypos,xpos]+=1
    n_observations[n_observations==0]=1
    acc_snowdepth/=n_observations
    
    return Bunch(snow=acc_snowdepth,x=xmap,y=ymap)


def write_data(data,filename):
    io.write(filename+"snowdepth",data.snow)
    io.write(filename+"x",data.x)
    io.write(filename+"y",data.y)

def main():
    files=glob.glob("??_2003")
    for f in files:
        os.chdir(f)
        print(f)
        xyzfiles=glob.glob("*BE*.xyz")
        for xf in xyzfiles:
            print(xf)
            data=load_data(xf)
            write_data(data,"../"+xf)
            
        os.chdir("../")

if __name__ == '__main__':
    main()