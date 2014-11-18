#!/usr/bin/env python
import os, glob

import numpy as np
import mygis

filesearch="{}_*_*_*_r*i*p*_*.nc"

geo_limits=dict(lat=[15,60],lon=[200,300])

def subset_file(inputfile,outputfile,limits):
    """docstring for subset_file"""
    fin=mygis.Dataset(inputfile)
    
    raise ValueError("This program has not been written yet")
    

def main(varname="tos"):
    """docstring for main"""
    files=glob.glob(filesearch.format(varname))
    if not os.path.exists("subset"):
        os.mkdir("subset")
        
    for f in files:
        print(f)
        subset_file(f,"subset/"+f,geo_limits)

if __name__ == '__main__':
    main()