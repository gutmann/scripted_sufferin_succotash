#!/usr/bin/env python
import glob
import re
import numpy as np
import swim_io
import matplotlib.pyplot as plt

def read_geosubset(filename,geo):
    lat=swim_io.read_nc(filename,"lat").data
    lon=swim_io.read_nc(filename,"lon").data-360
    goodlat=np.where((lat>=geo[0]) & (lat<=geo[1]))[0]
    goodlon=np.where((lon>=geo[2]) & (lon<=geo[3]))[0]
    subset=[goodlat[0],goodlat[-1]+1,goodlon[0],goodlon[-1]+1]
    return subset

def load_interannual(files,varname,subset=[0,None,0,None],timesubset=[0,None],geosubset=None):
    means=[]
    if geosubset:
        subset=read_geosubset(files[0],geosubset)
    for i,f in enumerate(files):
        d=swim_io.read_nc(f,var=varname)
        means.append(d.data[timesubset[0]:timesubset[1],subset[0]:subset[1],subset[2]:subset[3]].mean(axis=0))
        # d.ncfile.close()
    if re.match("SAR.*",files[0]):
        newmeans=[]
        for i in range(0,len(means),12):
            newmeans.append(np.dstack(means[i:i+12]).mean(axis=2))
        means=newmeans
    return np.dstack(means)# .std(axis=2)

def visualize_data(data,filename,varname):
    print(data.shape)
    if varname=="pr":
        data*=365 #make it annual precip not daily
        vmin=40.0
        vmax=400.0
    else:
        # vmin=data[data>0].min()
        # vmax=data[data<100].max()
        vmin=0
        vmax=4.0
    swim_io.write(filename.split('.')[0],data)
    data=data.std(axis=2)
    plt.clf()
    plt.imshow(data,vmin=vmin,vmax=vmax,origin="lower")
    plt.title(filename.split('.')[0].replace("_"," "))
    plt.colorbar()
    plt.draw()
    plt.savefig(filename)
    
    

def main():
    # 
    stats=["CAe1","SDe1","SDmon","CAe","SDe","CAe0","SDe0","CA","SD","SARe0"]#,"SARe1"]#,
            # "CAcold","CAhot","SDcold","SDhot","CAdry","CAwet","SDdry","SDwet"]
    # stats=["CA","SD"]
    stats=["CAe1","SDe1","SDmon","CAe","SDe","CAe0","SDe0","CA","SD","SARe0","SARe1"]
    # stats=["SARe0","SARe1"]
    bases=["narr","ncep"]
    vars=["pr","tasmin","tasmax"]
    resolutions=["12km","6km"]
    # resolutions=["12km"]
    BCstatus=["BC",""]
    yearsearch="200[1-8]"
    geosubset=None
    geosubset=[34,43.9,245-360,260-360]
    if geosubset==None:
        region="CONUS"
    else:
        region="hw"
    #stats=[]
    for s in stats:
        for b in bases:
            for v in vars:
                for r in resolutions:
                    for bc in BCstatus:
                        print(" ".join([s,b,v,r,bc]))
                        files=glob.glob(s+"/"+b+"/"+v+"/"+bc+s[:2]+"*"+r+"*"+yearsearch+"*.nc")
                        if len(files)>1:
                            data=load_interannual(files,v,geosubset=geosubset)
                            visualize_data(data,"_".join(["interannual",region,s,b,v,r,bc])+".png",v)
    # obs=["maurer.125","uw.4km","uw.0625"]
    if len(stats)<4:
        return
        
    obs=["maurer.125","uw.0625"]
    for v in vars:
        for o in obs:
        #for r in resolutions:
        #    if r=="12km":
        #        base_dir="../obs/maurer.125/"
        #    else:
        #        base_dir="../obs/uw.4km/"
            base_dir="../obs/"+o+"/"
            files=glob.glob(base_dir+v+"/*"+yearsearch+"*.nc")
            print(" ".join(["obs",o,v]))
            if len(files)>1:
                data=load_interannual(files,v,geosubset=geosubset)
                visualize_data(data,"_".join(["interannual",region,"obs",v,o])+".png",v)
    
if __name__ == '__main__':
    main()
