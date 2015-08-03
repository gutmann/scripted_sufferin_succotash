#!/usr/bin/env python

import sys
import glob

import numpy as np
import matplotlib.pyplot as plt

import load_data

outputfile="hist_{}.png"

def read_data(directory, ensemble):
    """docstring for read_data"""
    files=glob.glob(directory+"/*{}.txt".format(ensemble))
    outputdata=None
    for f in files:
        try:
            if outputdata==None:
                outputdata=load_data.cols(f)
            else:
                newdata=load_data.cols(f)
                outputdata=np.concatenate([newdata,outputdata])
        except ValueError:
            pass # file probably had nothing in it. 
    return outputdata
    
def main():
    """docstring for main"""
    files=glob.glob("current/cdp_ebr64_epac_????.txt")
    files.sort()
    
    for f in files:
        # plt.clf()
        ens=f.split("_")[-1][:4]
        print(ens)
        
        current=read_data("current",ensemble=ens)
        future=read_data("pgw",ensemble=ens)
        
        if (current.size>10):
            cdf,x,p=plt.hist(current.flat,bins=20,range=(0,10))
            print("Current: {}".format(100.*cdf[-1]/float(current.size)))
        if (future.size>10):
            cdf,x,p=plt.hist(future.flat,bins=20,range=(0,10))
            print("Future: {}".format(100.*cdf[-1]/float(current.size)))
            
        # if (current.size>10) and (future.size>10):
        #     plt.hist([current,future],bins=10,range=(5,10),label=["Current","PGW"],normed=True)
        #     plt.savefig(outputfile.format(ens))

if __name__ == '__main__':
    
    # ensemble=sys.argv[1]
    
    main()