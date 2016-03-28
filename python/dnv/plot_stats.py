#!/usr/bin/env python
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from dnv import load_wrf_tracks

def load_data(filesearch):
    files=glob(filesearch)
    files.sort()
    data=[]
    for f in files:
        keys,d = load_wrf_tracks.load_data(f)
        data.append(d)
    return keys, np.concatenate(data,axis=0)
    
ranges={"pmin":(-8000,-1000),
        "radius":(0,300),
        "wmean":(24,36),
        "wmax":(15,70),
        }

years = [2004, 2005, 2007, 2008, "????", "200[4,8]"]

def main():
    for y in years:
        print("")
        print("--------------------------------")
        print(y)
        keys, pgw_data  = load_data( "pgw/tracks_subset_{}.txt".format(y))
        keys, ctrl_data = load_data("ctrl/tracks_subset_{}.txt".format(y))
        
        output_filename="histograms_{}_{}.png"
        for i in range(len(keys)):
            if keys[i] in ranges:
                plt.clf()
                plt.hist( pgw_data[:,i][(np.isfinite( pgw_data[:,i]))], range=ranges[keys[i]], alpha=0.5, label="PGW")
                plt.hist(ctrl_data[:,i][(np.isfinite(ctrl_data[:,i]))], range=ranges[keys[i]], alpha=0.5, label="Current")
                plt.legend()
                plt.title(keys[i])
                plt.savefig(output_filename.format(keys[i], y))
                
                print(keys[i])
                if (keys[i]=="pmin"):
                    print("  pgw = {:8.5}".format(np.mean( pgw_data[:,i][(np.isfinite( pgw_data[:,i]))])))
                    print("  cur = {:8.5}".format(np.mean(ctrl_data[:,i][(np.isfinite(ctrl_data[:,i]))])))
                else:
                    print("  pgw = {:8.4}".format(np.mean( pgw_data[:,i][( pgw_data[:,i]>0)&(np.isfinite( pgw_data[:,i]))])))
                    print("  cur = {:8.4}".format(np.mean(ctrl_data[:,i][(ctrl_data[:,i]>0)&(np.isfinite(ctrl_data[:,i]))])))


if __name__ == '__main__':
    main()
