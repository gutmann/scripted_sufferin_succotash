#!/usr/bin/env python
import load_data
import matplotlib.pyplot as plt
import numpy as np


titles=["p99", "wetfrac0", "wetfrac1", "wetfrac2", "p95", "f95", "f99", "mean","meanfrac"]
lims=[10,0.05,0.05,0.05,5,0.05,0.05,150,0.3]
indicies=[1,2,3,0,4,7,5,6,7]

def plot_dataset(d):
    for i in range(2,19,2):
        index=indicies[i/2-1]
        plt.subplot(3,3,i/2)
        if i>17:
            dif=(d[:,index*2+2]-d[:,index*2+3])/d[:,index*2+2]
            index=8
        else:
            dif=d[:,index*2+2]-d[:,index*2+3]
        good=np.where(np.abs(dif)<900)[0]
        # var=np.std(dif[good])
        var=lims[index]
        # plt.clf();
        for cur in good:
            r=dif[cur]/var
            b=0.0
            if r<0: 
                b=np.abs(r)
                r=0.0
                if b>1:b=1
            if r>1: r=1
            plt.plot(d[cur,1],d[cur,0],'o',color=(r,0,b))
        plt.title(titles[index]+" Obs-WRF range= +/-"+str(var)[:5])
        plt.draw()

def main():
    snotel=load_data.cols("snotel_stats.txt")
    ghcn=load_data.cols("ghcn_comparison.txt")

    fig=plt.figure(figsize=(20,15),dpi=50)
    plot_dataset(snotel)
    plt.savefig("snotel_plots.png")
    plt.close()
    
    plt.figure(figsize=(20,15),dpi=50)
    plot_dataset(ghcn)
    plt.savefig("ghcn_plots.png")
    plt.close()
    

if __name__ == '__main__':
    main()