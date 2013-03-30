#!/usr/bin/env python
from glob import glob
import numpy as np
import swim_io
import os

alldirs=glob("agg*output")
alldirs.sort()
alldirs.append("subset_extremes")
distcheck="gamma*"
return_year=["  2yr", " 10yr", " 50yr","100yr"]
nyear=range(4)
drivingmodel="ncep"
drivingmodel="narr"
if len(sys.argv)>1:
    drivingmodel=sys.argv[1]

for y in nyear:
    for curdir in alldirs:
        os.chdir(curdir)
        files6=glob("*"+drivingmodel+"-pr-BC6km_values*"+distcheck)
        files6.sort()
        obsfile6=glob("extremes_obs-uw.0625-pr_values*"+distcheck)[0]
        obs6=swim_io.read_nc(obsfile6,"data").data
        d6=[]
        for f in files6:
            d6.append(swim_io.read_nc(f,"data").data)

        obsmean=np.mean(obs6[:,:,y])
        print("-"*20)
        print(" 6km "+str(curdir)+" "+return_year[y])
        print("method :   bias   %bias    RMS   %RMS    Obs")
        for i in range(len(files6)):
            tmp=np.where((d6[i][:,:,y]<200)&(d6[i][:,:,y]>1))
            diffs=d6[i][:,:,y][tmp]-obs6[:,:,y][tmp]
            method=files6[i].split("_")[1].split("-")[0]
            bias=np.mean(diffs)
            rms=np.sqrt(np.mean(diffs**2))
            print(" {0:6} : {1:6.1f} {2:6.1f}% {3:6.1f} {4:6.1f}% {5:6.1f}"
                    .format(method, bias, 100*bias/obsmean, rms, 100*rms/obsmean, obsmean))

        files12=glob("*"+drivingmodel+"-pr-BC12km_values*"+distcheck)
        files12.sort()
        obsfile12=glob("extremes_obs-maurer.125-pr_values*"+distcheck)[0]
        obs12=swim_io.read_nc(obsfile12,"data").data
        d12=[]
  
        for f in files12:
            d12.append(swim_io.read_nc(f,"data").data)

        obsmean=np.mean(obs12[:,:,y])
        print("-"*20)
        print(" 12km "+str(curdir)+" "+return_year[y])
        print("method :   bias   %bias    RMS   %RMS    Obs")
    
        for i in range(len(files12)):
            tmp=np.where((d12[i][:,:,y]<200)&(d12[i][:,:,y]>1))
            diffs=d12[i][:,:,y][tmp]-obs12[:,:,y][tmp]
            method=files12[i].split("_")[1].split("-")[0]
            bias=np.mean(diffs)
            rms=np.sqrt(np.mean(diffs**2))
            print(" {0:6} : {1:6.1f} {2:6.1f}% {3:6.1f} {4:6.1f}% {5:6.1f}"
                    .format(method, bias, 100*bias/obsmean, rms, 100*rms/obsmean, obsmean))
        os.chdir("../")
