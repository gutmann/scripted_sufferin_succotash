#!/usr/bin/env python
from stat_down import myio as io
import matplotlib.pyplot as plt
import numpy as np

def events(ds):
    storm_threshold=0.5
    nday_pre=5
    nday_post=5
    maxstorm=5
    
    instorm=False
    quietspell=0
    stormstart=0
    storm_length=0
    for i in range(len(ds)):
        # print(i,quietspell,ds[i],storm_threshold,instorm,storm_length)
        # if it didn't snow
        if ds[i]<storm_threshold:
            # add one to the count of days without snow
            quietspell+=1
            
        # if the count is greater than the pre-storm threshold (and we're not in a storm)
        if quietspell>=nday_pre and not instorm:
            # if it snowed this day, then we are in a storm!
            if ds[i]>storm_threshold:
                # record the storm start and set its length to 0
                instorm=True
                stormstart=i
                storm_length=0
                
        if instorm:
            # if we are already in a storm and it is still snowing, 
            # if ds[i]>storm_threshold:
                # add one to the storm length
            storm_length+=1
            # if this storm has lasted "too" long, reset statistics and find a new storm
            if storm_length>maxstorm:
                storm_length=0
                instorm=False
                quietspell=0
                
        # if we are "in" a storm and have hit enough quiet days after the storm
        if quietspell>=nday_post and instorm:
            # return this storms start point and length
            yield (stormstart,storm_length-nday_post)
            # we are no longer "in" a storm
            instorm=False
            storm_length=0
        if ds[i]>storm_threshold:
            # if it did snow set the count to 0
            quietspell=0
                

def main():
    # id=io.read_nc("~/dailySnotelDataWIn50kmGpsQC.nc","SNOTEL.station.id").data
    s=io.read_nc("~/dailySnotelDataWIn50kmGpsQC.nc","SWE").data
    d=io.read_nc("~/dailySnotelDataWIn50kmGpsQC.nc","snow.depth").data
    ds=np.diff(s)
    # dd=np.diff(d)
    d=d[1:]
    s=s[1:]
    event_list=[]
    for e,l in events(ds):
        # plt.plot(d[e-5:e+l+5])
        # plt.plot(s[e-5:e+l+5]*3)
        # plt.draw()
        # j=raw_input()
        
        storm_depth=d[e+l]-d[e]
        storm_size=s[e+l]-s[e]
        density=storm_size/storm_depth
        start=d[e+l]
        storm_collapse=[start-endpoint for endpoint in d[e+l+1:e+l+6]]
        event_list.append([density,storm_size,storm_depth,storm_collapse[0],storm_collapse[2],storm_collapse[-1],s[e+l+5]-s[e+l],e,l])
        
    return [event_list,s,d]
        

if __name__ == '__main__':
    main()

# import snow_checker
# reload(snow_checker)
# events=snow_checker.main()
# data=np.array(events[0])
# tmp=np.where((data[:,1]>2)&(data[:,2]>2)&(data[:,5]>0)&(np.abs(data[:,6])<1))
# plt.clf();plt.plot(data[tmp[0],0],data[tmp[0],5]/data[tmp[0],2],'x')