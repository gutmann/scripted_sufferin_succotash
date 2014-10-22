import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
# from matplotlib.colors import LightSource
from myLightSource import LightSource

cdict = {'red':   [(0.0,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],
         'green': [(0.0,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],
         'blue':  [(0.0,  1.0, 1.0),
                   (1.0,  1.0, 1.0)]}
whitemap=mpl.colors.LinearSegmentedColormap('white',cdict)


def shade(Z,colordata=None,cmap=None,clim=None,dx=1.0,dy=None,alt=65,az=100):
    """docstring for shade"""
    if dy==None:dy=dx
    if cmap==None:cmap=cm.terrain
    
    ls=LightSource(azdeg=az, altdeg=alt)
    shading=ls.shade(Z,cmap=whitemap,deltax=dx,deltay=dy)
    rgb=shading.copy()
    
    if colordata!=None:
        img=plt.imshow(colordata,cmap=cmap)
        if clim!=None:
            plt.clim(clim)
        img_data=img.to_rgba(colordata)
        
        hsv_shade=mpl.colors.rgb_to_hsv(shading[:,:,:3])
        hsv_color=mpl.colors.rgb_to_hsv(img_data[:,:,:3])
        
        #take the intensity data from the hill shade and the colors from the colordata
        hsv_color[:,:,2]=hsv_shade[:,:,2]
        
        # convert back to RGB
        rgb[:,:,:3]=mpl.colors.hsv_to_rgb(hsv_color)
        #combine the alpha channel from the Z and the colordata (e.g. masked elements)
        swap_alpha=(img_data[:,:,3]>rgb[:,:,3])
        rgb[:,:,3][swap_alpha]=img_data[:,:,3][swap_alpha]
    
    return rgb