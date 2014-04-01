#!/usr/bin/env python
import numpy as np
import mygis
import glob


def process(outputdata,data,index,func):
    """Combine two datasets on axis0[index] using func to decide which"""
    easy=np.where((data[index,...]>0) & (outputdata[index,...]==0))
    if len(easy[0])>0:
        outputdata[index,easy[0],easy[1]]=data[index,easy[0],easy[1]]
    
    hard=np.where((data[index,...]>0) & (outputdata[index,...]>0))
    for i in range(len(hard[0])):
        outputdata[index,hard[0][i],hard[1][i]]= \
            func([outputdata[index,hard[0][i],hard[1][i]],
                        data[index,hard[0][i],hard[1][i]]])
    

def main():
    """Combine a series of LIDAR processed netcdf files
    
    Files should have one variable "data" which is (2,y,x)
    assumes data[0,...] is a minimum and data[1,...] is a maximum
    When combining files, if two files have data for the same gridcell
    Take the min of data[0,...] from each and the max of data[1,...]
    
    Output a new file
    """
    outputfilename="compiled_lidar.nc"
    files=glob.glob("ot*.nc")
    outputdata=mygis.read_nc(files[0]).data
    for f in files[1:]:
        print(f)
        data=mygis.read_nc(f).data
        process(outputdata,data,index=0,func=np.min)
        process(outputdata,data,index=1,func=np.max)
        
    mygis.write(outputfilename,outputdata)

    
    

if __name__ == '__main__':
    main()