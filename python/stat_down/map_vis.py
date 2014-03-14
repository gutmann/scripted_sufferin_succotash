import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def vis(data,geo=[25,52.7,-124.7,-67],title="",vmin=None,vmax=None,proj='merc',
        cmap=None,colorbar=True,latstep=5.0,lonstep=10.0):
    """Plot a map of data using the bounds in geo=[lower_lat,upper_lat,left_lon,right_lon]
    
    Optionally specify a map title, min and max color values, colormap, projection,
    whether or not to draw a color bar, and what spacing to use between lat and lon grid lines
    """
    if geo=="subset":
        geo=[35,43,-113,-101]
    if geo=="conus":
        geo=[25,52.7,-124.7,-67]
        
    
    m = Basemap(projection=proj,llcrnrlat=geo[0],urcrnrlat=geo[1],\
                llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
    
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


