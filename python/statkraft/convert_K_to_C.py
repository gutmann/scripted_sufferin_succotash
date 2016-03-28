#!/usr/bin/env python

from netCDF4 import Dataset
import glob

def main():
    files=glob.glob("*/icar_*2d.nc")
    files.sort()
    
    for f in files:
        print(f)
        data=Dataset(f,mode='a')
        v=data.variables["ta2m"]
        v[:]=v[:]-273.15
        v.units="C"
        data.close()

if __name__ == '__main__':
    main()
