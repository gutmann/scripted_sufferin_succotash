#!/usr/bin/env python
import datetime,glob
import pickle
import numpy as np
from scipy.signal import medfilt2d
import matplotlib.pyplot as plt

import mygis as io
from bunch import Bunch

import sys
import flushprint
if flushprint.in_ipython():
    pass
else:
    sys.stdout=flushprint.Flushfile(sys.stdout)


def filter_dem(data):
    """Filter a LIDAR DEM to remove some spurious jumps"""
    return medfilt2d(data,5)

def calc_veg_height(veg_top,dem):
    veg1=veg_top-dem
    veg1[veg1<0]=0
    veg1[dem>3500]=0
    return medfilt2d(veg1,5)

def calc_snow_depth(snow_data,dem):
    return medfilt2d(snow_data,3)-dem


def decimate(data,factor=10, dfunc=np.mean):
    """Decimate data by a given factor in both dimension
    Uses averaging to decimate by default, could use median, min, max, ...
    """
    if factor<=1:
        return data

    sz=data.shape
    oned=(len(sz)==1)
    twod=(len(sz)==2)
    if twod:
        newsz=(np.floor(data.shape[0]/factor),np.floor(data.shape[1]/factor))
    elif oned:
        newsz=(np.floor(data.shape[0]/factor),)

    # the fast method if dfunc support axis= (e.g. most numpy stats: mean,min,max,median,...)
    try:
        if twod:
            data=data[:newsz[0]*factor,:newsz[1]*factor]
            data=data.reshape((newsz[0],factor,newsz[1],factor))
            outputdata=dfunc(dfunc(data,axis=3),axis=1)
        elif oned:
            data=data[:newsz[0]*factor]
            data=data.reshape((newsz[0],factor))
            outputdata=dfunc(data,axis=1)

    except TypeError:
    # slower but more general way incase dfunc doesn't support axis=  (should throw a TypeError)
        # newsz=(np.round(data.shape[0]/factor),np.round(data.shape[1]/factor))
        outputdata=np.zeros(newsz)
        if twod:
            for i in range(newsz[0]):
                for j in range(newsz[1]):
                    # assumes newsz*factor <= len(data)
                    outputdata[i,j]=dfunc(data[i*factor:(i+1)*factor,j*factor:(j+1)*factor])
        if oned:
            for i in range(newsz[0]):
                # assumes newsz*factor <= len(data)
                outputdata[i]=dfunc(data[i*factor:(i+1)*factor])

                # assumes len(data) < newsz*factor use this if newsz is computed with round instead of floor
                # if (i<(newsz[0]-1)):
                #     if (j<(newsz[1]-1)):
                #         outputdata[i,j]=dfunc(data[i*factor:(i+1)*factor,j*factor:(j+1)*factor])
                #     else:
                #         outputdata[i,j]=dfunc(data[i*factor:(i+1)*factor,j*factor:])
                # else:
                #     if (j<(newsz[1]-1)):
                #         outputdata[i,j]=dfunc(data[i*factor:,j*factor:(j+1)*factor])
                #     else:
                #         outputdata[i,j]=dfunc(data[i*factor:,j*factor:])

    return outputdata


def find_veg(vegz):
    """
    Find all locations within a grid of vegetation heights (vegz) that are primarily vegetation
    Decimates by a factor of 10x10 to get a 10m grid
    """
    return vegz>2

def load_lidar_geo(filename,decimation_factor=1):
    import pyproj
    ncdata=io.Dataset(filename)
    pvar=ncdata.variables["transverse_mercator"]
    # proj=str(pvar.spatial_ref)

    n=pvar.Northernmost_Northing
    s=pvar.Southernmost_Northing
    e=pvar.Easternmost_Easting
    w=pvar.Westernmost_Easting

    shape=ncdata.variables["Band1"].shape
    ncdata.close()
    dx=1.0
    xpoints=np.arange(w,e,dx)
    ypoints=np.arange(n,s,-dx)

    xpts=decimate(xpoints,decimation_factor)
    ypts=decimate(ypoints,decimation_factor)
    x,y=np.meshgrid(xpts,ypts)
    nx=xpts.shape[0]
    ny=ypts.shape[0]

    proj = pyproj.Proj(proj="utm",zone=13,ellps="WGS84")
    lon,lat = proj(x,y, inverse=True)
    # lat,lon=io.proj2ll(points=np.concatenate([x.reshape((-1,1)),y.reshape((-1,1))],axis=1),proj=proj)
    lat=lat.reshape((ny,nx))
    lon=lon.reshape((ny,nx))
    return Bunch(lat=lat,lon=lon)

def lidar2modscag(modscag_file="MODSCAG/fsca2008.dat",lidar_loc="lidar/",lidar_geo="snow-on-dem.nc",decimation_factor=5):
    from wsc import modscag,regrid

    print("loading modscag geo data")
    modscagdata=modscag.load(modscag_file)
    modscagdata.lat=modscagdata.lat[::-1]
    print("loading lidar data")
    data=load_fast(loc=lidar_loc,geofile=lidar_geo,decimation_factor=decimation_factor)
    if glob.glob("lidar2modscag.geolut.pickle"):
        print("reading geoLUT")
        geoLUT_depickler=pickle.Unpickler(open("lidar2modscag.geolut.pickle","r"))
        geoLUT=geoLUT_depickler.load()
    else:
        geoLUT=None

    print("computing mean")
    geoLUT,low_res_lidar_mean=regrid.agg(data,modscagdata.lat,modscagdata.lon,geo_lut=geoLUT,agg_func=np.mean)
    print("computing std")
    geoLUT,low_res_lidar_var =regrid.agg(data,modscagdata.lat,modscagdata.lon,geo_lut=geoLUT,agg_func=np.std)
    print("computing max")
    geoLUT,low_res_lidar_max =regrid.agg(data,modscagdata.lat,modscagdata.lon,geo_lut=geoLUT,agg_func=np.max)

    if not glob.glob("lidar2modscag.geolut.pickle"):
        print("writing geolut")
        geoLUT_pickler=pickle.Pickler(open("lidar2modscag.geolut.pickle","w"))
        geoLUT_pickler.dump(geoLUT)

    io.write("lidar2modscag_mean.nc", low_res_lidar_mean.data)
    io.write("lidar2modscag_std.nc", low_res_lidar_var.data)
    io.write("lidar2modscag_max.nc", low_res_lidar_max.data)




def load(snowoff="compiled_lidar.nc",snowon="snow-on-dem.nc"):
    """Load data from two netCDF files,
    convert data to snow depth, elevation, and veg height"""
    bare_data=io.read_nc(snowoff,"data").data
    dem=filter_dem(bare_data[0,...])
    veg=calc_veg_height(bare_data[1,...],dem)
    snow_data=io.read_nc(snowon,"Band1").data
    snow=calc_snow_depth(snow_data,dem)

    return Bunch(snow=snow,dem=dem,veg=veg,raw=[bare_data,snow_data])

def load_fast(dem="demf",snow="snowdepth",veg="vegmask",loc="./",geofile=None,decimation_factor=1):
    """docstring for load_fast"""
    demd=decimate(io.read_nc(loc+dem).data,decimation_factor)
    snowd=decimate(io.read_nc(loc+snow).data,decimation_factor)
    vegd=decimate(io.read_nc(loc+veg).data,decimation_factor)
    if geofile!=None:
        geo=load_lidar_geo(loc+geofile,decimation_factor=decimation_factor)
        lat=geo.lat
        lon=geo.lon
    else:
        print("No geofile supplied, try: snow-on-dem.nc?")
        lat=None
        lon=None
    lc=vegd.copy()
    lc[vegd<0.5]=2 #exposed
    lc[vegd>=0.5]=1#forest
    lc[demd<100]=0 #unclassified
    return Bunch(veg=np.round(vegd),lc=lc,data=snowd,snow=snowd, dem=demd,lat=lat,lon=lon,dates=[datetime.datetime(2010,5,5)])

# def bin_by_elevation(data,dem,mask,dz=50):
#     """docstring for bin_by_elevation"""
#     minz=dem[dem>100].min()
#     maxz=dem[dem<5000].max()
#
#     n=np.round((maxz-minz)/dz)
#     veg=np.zeros(n)
#     vegmed=np.zeros(n)
#     vegmin=np.zeros(n)
#     vegmax=np.zeros(n)
#     exposed=np.zeros(n)
#     exposedmed=np.zeros(n)
#     exposedmin=np.zeros(n)
#     exposedmax=np.zeros(n)
#     z=np.arange(minz,maxz+dz,dz)
#
#     for i in np.arange(n):
#         curexp=np.where((dem>z[i])&(dem<=z[i+1])&(mask==0)&np.isfinite(data)&(data>0)&(data<20))
#         curn=len(curexp[0])
#         if curn>0:
#             exposed[i]=np.mean(data[curexp])
#             sorted_data=np.sort(data[curexp])
#             exposedmin[i]=sorted_data[int(curn*0.1)]
#             exposedmed[i]=sorted_data[int(curn*0.5)]
#             exposedmax[i]=sorted_data[int(curn*0.9)]
#
#         curveg=np.where((dem>z[i])&(dem<=z[i+1])&(mask==1)&np.isfinite(data)&(data>0)&(data<20))
#         curn=len(curveg[0])
#         if curn>0:
#             veg[i]=np.mean(data[curveg])
#             sorted_data=np.sort(data[curveg])
#             vegmin[i]=sorted_data[int(curn*0.1)]
#             vegmed[i]=sorted_data[int(curn*0.5)]
#             vegmax[i]=sorted_data[int(curn*0.9)]
#
#     veg=np.ma.array(veg,mask=(veg==0))
#     vegmin=np.ma.array(vegmin,mask=(vegmin==0))
#     vegmed=np.ma.array(vegmed,mask=(vegmed==0))
#     vegmax=np.ma.array(vegmax,mask=(vegmax==0))
#
#     return Bunch(z=z[:n],veg=veg,vegmed=vegmed,vegmin=vegmin,vegmax=vegmax,
#                 exposed=exposed,exposedmed=exposedmed,exposedmin=exposedmin,exposedmax=exposedmax)


def main():
    import wsc.compare2lidar as c2l
    print("Loading data...")
    data=load_fast()

    decimation_factor=924 # 924 # for SNODAS or other 30arcsec grid
    decimation_factor=10
    print("Decimating Snow")
    snow=decimate(data.snow,decimation_factor)
    print("Decimating Veg")
    veg=decimate(data.veg,decimation_factor,dfunc=np.median)
    # veg_locations=find_veg(veg)
    print("Decimating DEM")
    dem=decimate(data.dem,decimation_factor)

    print("Binning")
    banded=c2l.bin_by_elevation(snow,dem,veg,dz=100)

    c2l.plot_elevation_bands(banded,"lidar_by_elev.png",title="Lidar snow depth at Niwot",
                        ylim=(0,6),ylabel="Snow Depth[m]")

    # print("Plotting")
    # plt.plot(banded.z,banded.exposed,label="Exposed",color="b",linewidth=2)
    # plt.plot(banded.z,banded.exposedmed,"--",label="Exp. Median",color="b",linewidth=2)
    # plt.bar([0],[1],color="skyblue",edgecolor="black",label="Exp. 10-90%")
    # plt.plot(banded.z,banded.veg,label="Vegetation",color="g",linewidth=2)
    # plt.plot(banded.z,banded.vegmed,"--",label="Veg. Median",color="g",linewidth=2)
    # plt.bar([0],[1],color="lightgreen",edgecolor="black",label="Veg. 10-90%")
    #
    # plt.plot(banded.z,banded.vegmin,color="black")
    # plt.plot(banded.z,banded.vegmax,color="black")
    # plt.plot(banded.z,banded.exposedmin,color="black")
    # plt.plot(banded.z,banded.exposedmax,color="black")
    #
    # plt.fill_between(banded.z,banded.exposedmin,banded.exposedmax,
    #                     color="skyblue",edgecolor="black")
    # plt.fill_between(banded.z,banded.vegmin,banded.vegmax,
    #                     color="lightgreen",alpha=0.5,edgecolor="black")
    #
    #
    # plt.legend(loc=2)
    # plt.xlim(2500,3800)
    # plt.ylim(0,6)
    # plt.ylabel("Snow Depth [m]")
    # plt.xlabel("Elevation [m]")
    # plt.title("Lidar snow depth at Niwot")
    # plt.savefig("lidar_by_elev.png")





if __name__ == '__main__':
    main()
