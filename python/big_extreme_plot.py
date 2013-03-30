#!/usr/bin/env python
import glob
import re

import numpy as np
import matplotlib.pyplot as plt

def get_from_file(filename,method,aggregation,resolution,nyear):
    with open(filename,"r") as f:
        l=""
        while not re.match(".*"+resolution+" *"+aggregation+" *"+nyear,l):
            l=f.next()
        while not re.match(".*"+method+" ",l):
            l=f.next()
    
    bias=float(l.split()[3][:-1])
    rms=float(l.split()[5][:-1])
    return (bias,rms)


def main():
    
    methods=["CA","CAe0","SDe0","SDmon","SARe0"]
    resolutions=["6km","12km"]
    aggs=["agg2output","agg4output","agg8output","subset_extremes"]
    aggs.reverse()
    year_intervals=["2yr","10yr","50yr","100yr"]
    bias=np.zeros((10,16))
    rms=np.zeros((10,16))
    filename="extremes_analysis.txt"
    
    i=0
    for y in year_intervals:
        for agg in aggs:
            j=0
            for m in methods:
                for res in resolutions:
                    b,r=get_from_file(filename,m,agg,res,y)
                    bias[j,i]=b
                    rms[j,i]=r
                    j+=1
            i+=1
    
    output=np.vstack([bias,rms])
    plt.imshow(output,vmin=-30,vmax=30)
    plt.title("Combined")
    plt.colorbar()
    plt.savefig("extreme_bigfig.png")

    plt.clf()
    plt.imshow(bias,vmin=-30,vmax=30)
    plt.title("Bias")
    plt.colorbar()
    plt.savefig("extreme_bias.png")

    plt.clf()
    plt.imshow(rms,vmin=0,vmax=50)
    plt.title("RMS error")
    plt.colorbar()
    plt.savefig("extreme_rms.png")

    
if __name__ == '__main__':
    main()
        