#!/usr/bin/env python
import os
import glob

import numpy as np
from stat_down import myio as io

def main():
    # test_files=["BCCA12K",  "BCSD12K",  "BCSDdisag12K",  "MAURER12K"]
    bd="/glade/scratch/gutmann/usbr/SW_tests/"
    test_files=glob.glob("*")
    for d in test_files:
        print(d)
        files=glob.glob(d+"/*Forcing*.nc")
        means=None
        n=0
        if len(files)>1:
            for f in files:
                try:
                    data=io.read_nc(f,"sw").data
                except Exception as e:
                    print(e,f)
                if means==None:
                    means=np.zeros(data.shape[1:])
                means+=data.mean(axis=0)
                n+=1
            if means!=None:
                means/=n
                print(bd+d+"_sw_mean")
                io.write(bd+d+"_sw_mean",means)

if __name__ == '__main__':
    main()