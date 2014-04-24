import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def vis(data,geo=[25,52.7,-124.7,-67],title="",vmin=None,vmax=None,proj='cyl',
        cmap=None,colorbar=True,latstep=5.0,lonstep=10.0,m=None):
    """Plot a map of data using the bounds in geo=[lower_lat,upper_lat,left_lon,right_lon]
    
    Optionally specify a map title, min and max color values, colormap, projection,
    whether or not to draw a color bar, and what spacing to use between lat and lon grid lines
    """
    if geo=="subset":
        if proj=="lcc":
            geo=[25.125,52.875,-124.75,-67]
        else:
            geo=[35,43,-113,-101]
    if geo=="conus":
        geo=[25,52.7,-124.7,-67]

    if not m:
        if proj=="cyl":
            m = Basemap(projection=proj,llcrnrlat=geo[0],urcrnrlat=geo[1],\
                        llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
        elif proj=="lcc":
            nx=data.shape[1]
            ny=data.shape[0]
            dx=4000.0
            m = Basemap(width=nx*dx,height=ny*dx,
                        rsphere=(6378137.00,6356752.3142),\
                        resolution='i',area_thresh=10000.,projection='lcc',\
                        lat_1=33.,lat_2=45.,lat_0=39.,lon_0=-107.0)
    
    
    mapimg=m.imshow(data,vmin=vmin,vmax=vmax,cmap=cmap)
    m.drawparallels(np.arange(20,60,latstep),labels=[1,0,0,0],dashes=[1,4])
    m.drawmeridians(np.arange(-120,-65,lonstep),labels=[0,0,0,1],dashes=[1,4])
    m.drawstates(linewidth=1.5)
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    
    if colorbar:
        m.colorbar()
    if title:
        plt.title(title)


