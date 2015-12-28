#!/usr/bin/env python
# [130:145,200:208]

import numpy as np
import mygis
import matplotlib.pyplot as plt

current_precc_file="b.e11.B20TRC5CNBDRD.f09_g16.{:03}.cam.h0.PRECC.192001-200512.nc"
current_precl_file="b.e11.B20TRC5CNBDRD.f09_g16.{:03}.cam.h0.PRECL.192001-200512.nc"
future_precc_file="b.e11.BRCP85C5CNBDRD.f09_g16.{:03}.cam.h0.PRECC.200601-208012.nc"
future_precl_file="b.e11.BRCP85C5CNBDRD.f09_g16.{:03}.cam.h0.PRECL.200601-208012.nc"

current_start_year=1920
future_start_year=2006

current_start_point = (1990-current_start_year)*12
current_end_point = (2000-current_start_year)*12

future_start_point = (2025-future_start_year)*12
future_end_point = (2035-future_start_year)*12


def main():
    """docstring for main"""
    plt.figure(figsize=(25,25))
    months=np.arange(1,13)
    for i in range(1,30):
        print(i)
        currentdata=(mygis.read_nc(current_precc_file.format(i+1),"PRECC",returnNCvar=True).data[current_start_point:current_end_point,130:145,200:208].mean(axis=1).mean(axis=1)
                    +mygis.read_nc(current_precl_file.format(i+1),"PRECL",returnNCvar=True).data[current_start_point:current_end_point,130:145,200:208].mean(axis=1).mean(axis=1))
        futuredata =(mygis.read_nc(future_precc_file.format(i+1), "PRECC",returnNCvar=True).data[future_start_point:future_end_point,  130:145,200:208].mean(axis=1).mean(axis=1)
                    +mygis.read_nc(future_precl_file.format(i+1), "PRECL",returnNCvar=True).data[future_start_point:future_end_point,  130:145,200:208].mean(axis=1).mean(axis=1))
        
        cur=np.reshape(currentdata,(10,12)).mean(axis=0)
        fut=np.reshape(futuredata,(10,12)).mean(axis=0)
        delta=1000.*86400.*(fut-cur)
        
        plt.subplot(6,5,i+1)
        # plt.plot(months,delta,label="{:03}".format(i+1),color="black")
        for j in range(12):
            if delta[j]>0:
                plt.fill_between([months[j]-0.4,months[j]+0.4],[delta[j],delta[j]],color="blue")
            else:
                plt.fill_between([months[j]-0.4,months[j]+0.4],[delta[j],delta[j]],color="red")
        plt.plot([0,13],[0,0],":",color="black")
        plt.ylim(-0.5,0.5)
        plt.xlim(0.1,12.9)
        plt.title("Ens {:03}".format(i+1))
        
    # plt.legend()
    plt.tight_layout
    plt.savefig("all_diffs.png")
        

if __name__ == '__main__':
    main()