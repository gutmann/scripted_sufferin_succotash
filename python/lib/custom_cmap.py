import matplotlib.pyplot as plt
import numpy as np

import matplotlib.colors as mcolors


def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    
    from : http://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
    
    example usage: 
    c = mcolors.ColorConverter().to_rgb
    rvb = make_colormap(
        [c('red'), c('violet'), 0.33, c('violet'), c('blue'), 0.66, c('blue')])
    
    """
    # this might make it use black at the bottom and white at the very top, though it is transitioning at 0 and 1 so...
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            if len(seq[i - 1])==3:
                r1, g1, b1 = seq[i - 1]
            else:
                r1, g1, b1, a1 = seq[i - 1]
            if len(seq[i + 1])==3:
                r2, g2, b2 = seq[i + 1]
            else:
                r2, g2, b2, a2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


def subset(cmap=None,clim=(0,255),step=1):
    """docstring for subset"""
    if type(cmap)==str:
        cmap=plt.cm.__getattribute__(cmap)
    
    newcmap=[]
    for i in range(clim[0],clim[1]+step,step):
        newcmap.append(cmap(i))
        if (i<clim[1]):
            newcmap.append(float(i-clim[0])/(clim[1]-clim[0]))
            
    return make_colormap(newcmap)
    
    
    
    