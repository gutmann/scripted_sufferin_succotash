#!/usr/bin/env python
import glob,os,sys
import os.path
import numpy as np
# import swim_io as io
from bunch import Bunch
import mygis

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
    mygis.write(filename+"snowdepth",data.snow)
    mygis.write(filename+"x",data.x)
    mygis.write(filename+"y",data.y)

def clpx():
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

def niwot_grid():
    upper_left=np.array([438999.5,4436000.5])
    dx=1.0
    dy=-1.0
    nx=21001
    ny=18001
    x=(upper_left[0]+dx/2) + np.arange(nx)*dx
    y=(upper_left[1]+dy/2) + np.arange(ny)*dy
    return x,y
    
def slow_process(xyzfile,data,x,y):
    with open(xyzfile,'ru') as f:
        for l in f:
            curline=l.split()
            curx=np.float(curline[0])
            cury=np.float(curline[1])
            curz=np.float(curline[2])
            xpos=np.argmin(np.abs(x-curx))
            ypos=np.argmin(np.abs(y-cury))
            if data[0,ypos,xpos]==0:
                data[:,ypos,xpos]=curz
            else:
                if curz<data[0,ypos,xpos]:
                    data[0,ypos,xpos]=curz
                elif curz>data[1,ypos,xpos]:
                    data[1,ypos,xpos]=curz

def fast_process(xyzfile,data,x,y):
    filedata=np.loadtxt(xyzfile)
    for curline in filedata:
        curx=curline[0]
        cury=curline[1]
        curz=curline[2]
        xpos=np.argmin(np.abs(x-curx))
        ypos=np.argmin(np.abs(y-cury))
        if data[0,ypos,xpos]==0:
            data[:,ypos,xpos]=curz
        else:
            if curz<data[0,ypos,xpos]:
                data[0,ypos,xpos]=curz
            elif curz>data[1,ypos,xpos]:
                data[1,ypos,xpos]=curz


def niwot():
    # xyzfile="giant_text_file.txt"
    # xyzfile="ot_439000_4418000.txt"
    xyzfiles=sys.argv[1:]

    demfile=xyzfiles[0].replace(".txt",".nc")
    print(demfile)
    if os.path.isfile(demfile):
        print("Already found:"+demfile)
        return
    else:
        mygis.write(demfile,np.zeros((10,10)))
    
    x,y=niwot_grid()
    data=np.zeros((2,y.size,x.size),dtype="f")
    success=False
    for curfile in xyzfiles:
        print(curfile)
        try:
            fast_process(curfile,data,x,y)
            success=True
        except:
            print("Slower processing")
            try:
                slow_process(curfile,data,x,y)
                success=True
            except:
                print("All processing failed for file: "+curfile)
    
    # if os.path.isfile(demfile):
    #     os.remove(demfile)
    if success:
        try:
            os.remove(demfile)
        except:
            pass
        # should automatically clobber the temporary file anyway
        mygis.write(demfile,data)
    

if __name__ == '__main__':
    # clpx()
    print("Running for Niwot Lidar data")
    niwot()