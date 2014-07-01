import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def vis(data,geo=[25,52.7,-124.7,-67],title="",vmin=None,vmax=None,ylim=None,xlim=None,proj='cyl',
        cmap=None,colorbar=True,latstep=5.0,lonstep=10.0,m=None,lat=None,lon=None,clim=None,
        reproject=False,width=None,height=None):
    """Plot a map of data using the bounds in geo=[lower_lat,upper_lat,left_lon,right_lon]
    
    Optionally specify a map title, min and max color values, colormap, projection,
    whether or not to draw a color bar, and what spacing to use between lat and lon grid lines
    """
    if not m:
        if geo=="subset":
            if proj=="lcc":
                pass # we don't actually need geo for the lcc (WRF) subset
                # geo=[25.125,52.875,-124.75,-67]
            else:
                geo=[35,43,-113,-101]
        if geo=="conus":
            geo=[25,52.7,-124.7,-67]
    
        
        if geo==None:
            if (lat==None) or (lon==None):
                raise TypeError("Missing either geo or lat/lon variables")
        
            georank=len(lat.shape)
            if georank>1:
                if georank>2:
                    lat=lat[0,...]
                    lon=lon[0,...]
                dlat=lat[1,0]-lat[0,0]
                dlon=lon[0,1]-lon[0,0]
                lat0=lat[0,0]
                lon0=lon[0,0]
                lat1=lat[-1,0]
                lon1=lon[0,-1]
            else:
                dlat=lat[1]-lat[0]
                dlon=lon[1]-lon[0]
                lat0=lat[0,0]-dlat/2
                lon0=lon[0,0]-dlon/2
                lat1=lat[-1,0]+dlat/2
                lon1=lon[0,-1]+dlon/2
            geo=[lat0,lat1,lon0,lon1]

        if proj=="cyl":
            m = Basemap(projection=proj,llcrnrlat=geo[0],urcrnrlat=geo[1],\
                        llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
        elif proj=="lcc":
            nx=data.shape[1]
            ny=data.shape[0]
            dx=4000.0
            if not width:
                width=nx*dx
                height=ny*dx
            m = Basemap(width=width,height=height,
                        rsphere=(6378137.00,6356752.3142),\
                        resolution='i',area_thresh=10000.,projection='lcc',\
                        lat_1=33.,lat_2=45.,lat_0=39.,lon_0=-107.0)


    if clim:
        vmin=clim[0]
        vmax=clim[1]
    
    if reproject:
        #note, to reproject lat and lon must be specified
        if len(lon.shape)==1:
            x, y = m(*np.meshgrid(lon,lat))
        else:
            x, y = m(lon,lat)
        cs = m.pcolormesh(x,y,data,vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        mapimg=m.imshow(data,vmin=vmin,vmax=vmax,cmap=cmap)
    m.drawparallels(np.arange(20,60,latstep),labels=[1,0,0,0],dashes=[1,4])
    m.drawmeridians(np.arange(-120,-65,lonstep),labels=[0,0,0,1],dashes=[1,4])
    m.drawstates(linewidth=1.5)
    m.drawcountries(linewidth=1.5)
    m.drawcoastlines(linewidth=1.5)
    
    if colorbar:
        m.colorbar()
    if title:
        plt.title(title)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    
    return m


