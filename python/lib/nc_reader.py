#!/usr/bin/env python

"""
SYNOPSIS

    NC_Reader class makes it "easy" to read input from one netCDF file and match it to another
    
DESCRIPTION

    Library routines for reading a netCDF dataset and interpolating
    
EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION
    0.1
    
"""

import numpy as np

import swim_io
from bunch import Bunch
import glob

# given a set of 4 points (x[0-3],y[0-3]) and an internal point (xi,yi)
#  calculate the weights to assigne to each of the four points to bilinearly
#  interpolate their z values to the point xi,yi
def bilin_weights(yi,y,xi,x):
    x0=np.abs((xi-x[0])/(x[1]-x[0]))
    x1=1-x0 # equivalent to: np.abs((xi-x[1])/(x[1]-x[0]))
    x2=np.abs((xi-x[2])/(x[3]-x[2]))
    x3=1-x2 # equivalent to: np.abs((xi-x[3])/(x[3]-x[2]))
    y5=y[0]*x1+y[1]*x0
    y6=y[2]*x3+y[3]*x2
    f1=(yi-y5)/(y6-y5)
    f2=1-f1 # equivalent to: (y6-yi)/(y6-y5)
    return np.array([x1*f2,x0*f2,x3*f1,x2*f1])


def match_xy(lat1,lon1,lat2,lon2,subset=None,scheme="bilin"):
    '''Create a Geographic Look Up Table
    
    lat/lon inputs should be grids of latitude and longitude
        lat1/lon1 = low resolution input grid
        lat2/lon2 = high resolution output grid
    '''
    # (lon2,lat2)=read_geo_latlon(georef,subset=subset)
    # if lon1.max()>180:
    #     lon1=lon1-360
    # if lon2.max()>180:
    #     lon2=lon2-360
    N=lat2.shape
    N1=lat1.shape
    # output data
    geoLUT=np.empty((N[0],N[1],4,3))

    # intermediate variables
    x=np.zeros(4).astype('i')
    y=np.zeros(4).astype('i')
    winhalfsize=3 # search window for next data point
    
    # figure out which direction is positive in latitude and longitude (just in case)
    dxinc=np.sign(lon1[1,1]-lon1[0,0]).astype('i')
    dyinc=np.sign(lat1[1,1]-lat1[0,0]).astype('i')

    # loop over y
    for i in range(N[1]):
        # find the first positions so the rest can be done faster relative to it
        j=0
        # compute distances for the ENTIRE low-res grid for the first input point
        dists=(lat1-lat2[j,i])**2 + (lon1-lon2[j,i])**2
        # find the minimum distance position
        (lasty,lastx)=np.unravel_index(dists.argmin(), dists.shape)
        (prevx,prevy)=(lastx,lasty)
        # create a window around that position to start searching as we loop through the rest of it
        y0=max(lasty-winhalfsize,0)
        y1=min(lasty+winhalfsize,N1[0])
        x0=max(lastx-winhalfsize,0)
        x1=min(lastx+winhalfsize,N1[1])
        latwin=lat1[y0:y1,x0:x1]
        lonwin=lon1[y0:y1,x0:x1]
        for j in range(N[0]):
            # if we have moved update the window position
            if (prevx!=lastx) or (prevy!=lasty):
                y0=max(lasty-winhalfsize,0)
                y1=min(lasty+winhalfsize,N1[0])
                x0=max(lastx-winhalfsize,0)
                x1=min(lastx+winhalfsize,N1[1])
                latwin=lat1[y0:y1,x0:x1]
                lonwin=lon1[y0:y1,x0:x1]
                (prevx,prevy)=(lastx,lasty)
            # compute distances just over this window
            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
            # find minimum distance position
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            # convert that position to a position in the large grid
            lastx=newx+x0
            lasty=newy+y0
            # if i>190:
            #     print(newx,newy)
            #     print(x0,x1,y0,y1)
            #     print(lastx,lasty)
            # this is the first point we care about
            x[0]=newx
            y[0]=newy
            # now find the position of the other three surrounding points
            if scheme!="nn":
                winsz=latwin.shape
                if latwin[newy,newx]<lat2[j,i]:
                    if lonwin[newy,newx]<lon2[j,i]:
                        x[1]=max(min(newx+dxinc,winsz[1]-1),0)
                        x[2]=newx
                        x[3]=max(min(newx+dxinc,winsz[1]-1),0)
                        y[1]=newy
                        y[2]=max(min(newy+dyinc,winsz[0]-1),0)
                        y[3]=max(min(newy+dyinc,winsz[0]-1),0)
                    else:
                        x[1]=max(min(newx-dxinc,winsz[1]-1),0)
                        x[2]=newx
                        x[3]=max(min(newx-dxinc,winsz[1]-1),0)
                        y[1]=newy
                        y[2]=max(min(newy+dyinc,winsz[0]-1),0)
                        y[3]=max(min(newy+dyinc,winsz[0]-1),0)
                else:
                    if lonwin[newy,newx]<lon2[j,i]:
                        x[1]=max(min(newx+dxinc,winsz[1]-1),0)
                        x[2]=newx
                        x[3]=max(min(newx+dxinc,winsz[1]-1),0)
                        y[1]=newy
                        y[2]=max(min(newy-dyinc,winsz[0]-1),0)
                        y[3]=max(min(newy-dyinc,winsz[0]-1),0)
                    else:
                        x[1]=max(min(newx-dxinc,winsz[1]-1),0)
                        x[2]=newx
                        x[3]=max(min(newx-dxinc,winsz[1]-1),0)
                        y[1]=newy
                        y[2]=max(min(newy-dyinc,winsz[0]-1),0)
                        y[3]=max(min(newy-dyinc,winsz[0]-1),0)
                # finally compute the weights for each of the four surrounding points for a bilinear interpolation
                # simple inverse distance weighting
                if scheme=="inverse":
                    dists=1.0/dists
                    weights=dists[y,x]/np.sum(dists[y,x])

                # bilinear interpolation if the grid is perfectly regular (most Atm. model grids are NOT in lat/lon space)
                if scheme=="regular":
                    denom=(np.abs(latwin[y[0],x[0]]-latwin[y[1],x[1]])*np.abs(lonwin[y[0],x[0]]-lonwin[y[1],x[1]]))
                    w1=((np.abs(latwin[y[1],x[1]]-lat2[j,i])*np.abs(lonwin[y[1],x[1]]-lon2[j,i]))/denom)
                    w2=((np.abs(latwin[y[0],x[0]]-lat2[j,i])*np.abs(lonwin[y[0],x[0]]-lon2[j,i]))/denom)
                    w3=((np.abs(latwin[y[0],x[1]]-lat2[j,i])*np.abs(lonwin[y[0],x[1]]-lon2[j,i]))/denom)
                    w4=((np.abs(latwin[y[1],x[0]]-lat2[j,i])*np.abs(lonwin[y[1],x[0]]-lon2[j,i]))/denom)
                    weights=[w1,w2,w3,w4]
                
                # bilinear interpolation for an arbitrary grid spacing 
                # (must be a grid, but this will handle a lot of irregularities)
                if scheme=="bilin":
                    weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x])
                
                # store the current results in the output array
                geoLUT[j,i,:,0]=y+y0
                geoLUT[j,i,:,1]=x+x0
                geoLUT[j,i,:,2]=weights
            else:
                geoLUT[j,i,0,0]=y[0]+y0
                geoLUT[j,i,0,1]=x[0]+x0
    
    if scheme=="nn":
        return (geoLUT[:,:,0,1],geoLUT[:,:,0,0])
    else:
        return geoLUT


class NC_Reader(object):
    _filenames=None
    _curfile=0
    x=-99
    y=-99
    _bilin=False
    _nn=True
    _geoLUT=None
    posinfile=0
    npos=1
    latvarname="latitude"
    lonvarname="longitude"
    _vars=list(["RAINNC"])
    xmin=0;xmax=None;ymin=0;ymax=None
    # xmin=100;xmax=-100;ymin=100;ymax=-100
    
    def init_xy(self, geoin_file,geomatch_file,subset=None,geo_subset=None,glatvar="XLAT",glonvar="XLONG"):
        # narrfilename='/Volumes/Data2/usbr/NARR/flxnc/HGT_2006010100.nc'
        # wrffilename='/Volumes/G-SAFE/headwaters/wrf_output/terrain/2km_wrf_input_d01'
        nlat=swim_io.read_nc(geoin_file,var=self.latvarname).data
        nlon=swim_io.read_nc(geoin_file,var=self.lonvarname).data# -360
        if nlon.max()>180:
            nlon-=360
        if len(nlon.shape)==1:
            nlon,nlat=np.meshgrid(nlon,nlat)
        latinfo=swim_io.read_nc(geomatch_file,var=glatvar,returnNCvar=True)
        if len(latinfo.data.shape)==3:
            wlat=latinfo.data[0,self.ymin:self.ymax,self.xmin:self.xmax]
        elif len(latinfo.data.shape)==2:
            wlat=latinfo.data[self.ymin:self.ymax,self.xmin:self.xmax]
        else: # len(latinfo.data.shape)==1:
            wlat=latinfo.data[self.ymin:self.ymax]
        latinfo.ncfile.close()
        loninfo=swim_io.read_nc(geomatch_file,var=glonvar,returnNCvar=True)
        if len(loninfo.data.shape)==3:
            wlon=loninfo.data[0,self.ymin:self.ymax,self.xmin:self.xmax]
        elif len(loninfo.data.shape)==2:
            wlon=loninfo.data[self.ymin:self.ymax,self.xmin:self.xmax]
        else: # len(loninfo.data.shape)==1:
            wlon=loninfo.data[self.xmin:self.xmax]
            wlon,wlat=np.meshgrid(wlon,wlat)
        if wlon.max()>180:
            wlon-=360
        loninfo.ncfile.close()
        if geo_subset:
            ymin=np.where(wlat>=geo_subset[0])[0][0]
            ymax=np.where(wlat<=geo_subset[1])[0][-1]
            xmin=np.where(wlon>=geo_subset[2])[1][0]
            xmax=np.where(wlon<=geo_subset[3])[1][-1]
            subset=[ymin,ymax,xmin,xmax]
        if subset:
            wlat=wlat[subset[0]:subset[1],subset[2]:subset[3]]
            wlon=wlon[subset[0]:subset[1],subset[2]:subset[3]]
        if self._nn:
            (x,y)=match_xy(nlat,nlon,wlat,wlon,scheme="nn")
            self.x=x.astype('i')
            self.y=y.astype('i')
        elif self._bilin:
            self._geoLUT=match_xy(nlat,nlon,wlat,wlon,scheme="bilin")
        
        
    def __init__(self, file_search,geoin_file=None,geomatch_file=None,ntimes=1,readvars=list(['RAINNC']),
                    nn=True, bilin=False, latvar=None,lonvar=None,firstfile_timeinit=0,subset=None,geo_subset=None,
                    glatvar="XLAT",glonvar="XLONG",filelist=None,*args, **kwargs):
        super(NC_Reader,self).__init__()# *args, **kwargs)
        if filelist==None:
            self._filenames=glob.glob(file_search)
        else:
            self._filenames=filelist
        self._filenames=np.sort(self._filenames)
        if latvar!=None:
            self.latvarname=latvar
        if lonvar!=None:
            self.lonvarname=lonvar
            
        if bilin:
            self._nn=False
            self._bilin=True
        else:
            self._nn=True
            self._bilin=False
        if ntimes>1:
            d=swim_io.read_nc(self._filenames[0],var=readvars[0],returnNCvar=True)
            ntimes=d.data.shape[0]
            d.ncfile.close()
            self.npos=ntimes
        self.posinfile=firstfile_timeinit
        if geomatch_file:
            if geoin_file==None:
                geoin_file=self._filenames[0]
            # print('Calculating XY match lookup table, this may take a while.')
            self.init_xy(geoin_file,geomatch_file,subset=subset,geo_subset=geo_subset,
                        glatvar=glatvar,glonvar=glonvar)
        else:
            d=swim_io.read_nc(self._filenames[0],var=readvars[0],returnNCvar=True)
            if geo_subset:
                lat=swim_io.read_nc(self._filenames[0],var=self.latvarname).data
                lon=swim_io.read_nc(self._filenames[0],var=self.lonvarname).data
                if lon.max()>180:
                    lon-=360
                if len(lat.shape)==1:
                    ymin=np.where(lat>=geo_subset[0])[0][0]
                    ymax=np.where(lat<=geo_subset[1])[0][-1]
                    xmin=np.where(lon>=geo_subset[2])[0][0]
                    xmax=np.where(lon<=geo_subset[3])[0][-1]
                else:
                    ymin=np.where(lat>=geo_subset[0])[0][0]
                    ymax=np.where(lat<=geo_subset[1])[0][-1]
                    xmin=np.where(lon>=geo_subset[2])[1][0]
                    xmax=np.where(lon<=geo_subset[3])[1][-1]
                subset=[ymin,ymax,xmin,xmax]

            shp=d.data.shape
            if len(shp)==2:
                (ny,nx)=d.data.shape
            else:
                (ntimes,ny,nx)=d.data.shape
            d.ncfile.close()
            if subset:
                y=np.arange(ny)[subset[0]:subset[1]]
                x=np.arange(nx)[subset[2]:subset[3]]
            else:
                y=np.arange(ny)
                x=np.arange(nx)
            (x,y)=np.meshgrid(x,y)
            self.x=x
            self.y=y
        self._vars=readvars
    
    
    # we are our own iterator...
    def __iter__(self):
        return self

    def nextNN(self):
        curfile=self._curfile
        if curfile>=len(self._filenames):
            raise StopIteration
        minx=self.x.min()
        maxx=self.x.max()+1
        miny=self.y.min()
        maxy=self.y.max()+1
        curx=self.x-minx
        cury=self.y-miny
        data=list()
        for thisvar in self._vars:
            ncfile=swim_io.read_nc(self._filenames[curfile],var=thisvar,returnNCvar=True)
            if self.npos==1:
                data.append(ncfile.data[miny:maxy,minx:maxx][cury,curx])
            else:
                data.append(ncfile.data[self.posinfile,miny:maxy,minx:maxx][cury,curx])
            ncfile.ncfile.close()
        self.posinfile+=1

        if (self.npos==1) or (self.posinfile>=self.npos):
            self._curfile+=1
            self.posinfile=0
            if self.npos>1:
                if self._curfile<len(self._filenames):
                    d=swim_io.read_nc(self._filenames[self._curfile],var=self._vars[0],returnNCvar=True)
                    ntimes=d.data.shape[0]
                    d.ncfile.close()
                    self.npos=ntimes
            
        
        return data
    

    def nextBilin(self):
        curfile=self._curfile
        if curfile>=len(self._filenames):
            raise StopIteration
        geoLUT=self._geoLUT
        x=geoLUT[:,:,:,1]
        y=geoLUT[:,:,:,0]
        w=geoLUT[:,:,:,2]
        minx=x.min().astype('i')
        maxx=x.max().astype('i')+1
        miny=y.min().astype('i')
        maxy=y.max().astype('i')+1
        curx=(x-minx).astype('i')
        cury=(y-miny).astype('i')
        N=curx[:,:,0].shape
        data=list()
        for thisvar in self._vars:
            thisdata=np.zeros(N)
            ncfile=swim_io.read_nc(self._filenames[curfile],var=thisvar,returnNCvar=True)
            if self.npos==1:
                curdata=ncfile.data[miny:maxy,minx:maxx]
            else:
                curdata=ncfile.data[self.posinfile,miny:maxy,minx:maxx]
                self.posinfile+=1
            for i in range(4):thisdata+=curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            data.append(thisdata)
            ncfile.ncfile.close()
        
        # datestr=self._filenames[curfile].split('.')[1]

        if (self.npos==1) or (self.posinfile>=self.npos):
            self._curfile+=1
            self.posinfile=0

        return data
        
    def next(self):
        if self._geoLUT==None:
            return self.nextNN()
        else:
            return self.nextBilin()
    # self.__next__=self.next
            
    def close(self):
        pass
    def __enter__(self):
        return self
    
    def __exit__(self):
        self.close()
    
    def __del__(self):
        self.close()
        
