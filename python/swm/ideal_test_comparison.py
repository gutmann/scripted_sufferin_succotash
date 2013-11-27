#!/usr/bin/env python
from glob import glob

import numpy as np
import matplotlib.pyplot as plt

from stat_down import myio as io



def main():
    """Compare data from various ideal swim tests"""

    directories=glob("output_*")
    for d in directories:
        files=glob(d+"/*")
        files.sort()
        print(d)
        if len(files)>30:
            rain1=io.read_nc(files[-3],"rain").data
            rain2=io.read_nc(files[-2],"rain").data
            ny,nx=rain1.shape
            rainrate=(rain2[ny/2,:]-rain1[ny/2,:])/3600.
            x=np.arange(nx)/nx*2
            plt.clf();
            plt.plot(rainrate)
            plt.ylim(0,0.001)
            plt.title(d)
            plt.xlabel("Distance [km]")
            plt.ylabel("Precipitation rate [mm/s]")
            plt.savefig(d+".png")
            
        
        
    

if __name__ == '__main__':
    main()
