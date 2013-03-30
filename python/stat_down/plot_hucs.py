#!/usr/bin/env python
import glob

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import shapefile

def plot2map(shape,m,mapbox,color="blue",thick=1.0,label=None):

    plotlabel=True
    for poly in shape:
        bbox=poly.shape.bbox
        parts=poly.shape.parts
        points=poly.shape.points
    
        # inbox=np.where((lon>bbox[0]) & (lon<=bbox[2]) & (lat>=bbox[1]) & (lat<=bbox[3]))
        # print((bbox[2],mapbox[2]) , (bbox[0],mapbox[3]) , (bbox[3],mapbox[0]) , (bbox[1],mapbox[1]))
        # print((bbox[2]<mapbox[2]) , (bbox[0]>mapbox[3]) , (bbox[3]<mapbox[0]) , (bbox[1]>mapbox[1]))
        if (bbox[2]<mapbox[2]) or (bbox[0]>mapbox[3]) or (bbox[3]<mapbox[0]) or (bbox[1]>mapbox[1]):
            pass
        else:
            for i in range(len(parts)):
                if i==len(parts)-1:
                    xypoints=np.array(points[parts[i]:])
                else:
                    xypoints=np.array(points[parts[i]:parts[i+1]])
                if plotlabel:
                    m.plot(xypoints[:,0],xypoints[:,1],color=color,linewidth=thick,label=label)
                    plotlabel=False
                else:
                    m.plot(xypoints[:,0],xypoints[:,1],color=color,linewidth=thick)


def main():
    files=glob.glob("HUC*/*.shp")
    files.sort()
    # files=[files[0]]
    files=files[:3]
    print(files)
    conus=True
    colors=["red","blue","green","magenta","cyan"]
    thicknesses=[2.5,2,1.5,1,0.5]
    geo=[35,43,-113,-101]
    if conus:
        geo=[25,52,-125,-65]
        thicknesses=[0.5,0.3,0.2]
        plt.figure(figsize=(10,5))
    m = Basemap(projection='cyl',llcrnrlat=geo[0],urcrnrlat=geo[1],\
                llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
    for i in range(len(files)):
        curpos=len(files)-i-1
        filename=files[curpos]
        print(filename)
        shp=shapefile.Reader(filename)
        polygons=shp.shapeRecords()
        plot2map(polygons,m,geo,color=colors[curpos],thick=thicknesses[curpos],label=filename.split("/")[0])
        
    if conus:
        m.drawparallels(np.arange(25,55,5.),labels=[1,0,0,0],dashes=[1,4])
        m.drawmeridians(np.arange(-120,-65,10.),labels=[0,0,0,1],dashes=[1,4])
        m.drawstates(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawcoastlines(linewidth=0.5)
    else:
        m.drawstates(linewidth=1.5,color='k')
        m.drawparallels(np.arange(36,43,2.),labels=[1,0,0,0],dashes=[1,4])
        m.drawmeridians(np.arange(-112,-103,4.),labels=[0,0,0,1],dashes=[1,4])
    plt.legend(loc=4)
    if conus:
        region="conus"
    else:
        region="subdomain"
    plt.savefig(region+"_hucmap.png",dpi=200)


if __name__ == '__main__':
    main()