import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def vis(data,geo=[],title="",vmin=None,vmax=None,cmap=None):
    """Plot a map of data using the bounds in geo=[lllat,urlat,lllon,urlon]
    
    Optionally specify a map title, min and max value and colormap
    """
    # geo=[35,43,-113,-101]
    m = Basemap(projection='cyl',llcrnrlat=geo[0],urcrnrlat=geo[1],\
                llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
    
    mapimg=m.imshow(data,vmin=vmin,vmax=vmax,cmap=cmap)
    # m.drawparallels(np.arange(25,55,5.),labels=[1,0,0,0],dashes=[1,4])
    # m.drawmeridians(np.arange(-120,-65,10.),labels=[0,0,0,1],dashes=[1,4])
    m.drawstates(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    
    m.colorbar()
    plt.title(title)
