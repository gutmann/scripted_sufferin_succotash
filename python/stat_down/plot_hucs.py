#!/usr/bin/env python
import glob

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import shapefile

huc_location="/d2/gutmann/usbr/hucdata/"

geosubsets=dict(
    cohw=[32.0625, 44.9375, -114.9375, -99.0625], #[35,43,-112.8,-101.7], # subdomain # 32.0625, 44.9375, -114.9375, -99.0625
    pnw=[39.0625, 50, -124.6875, -108.8125], #N is 52.8125
    # se=[25.1875, 39.0625, -93.6875, -75], #E is -67.0625
    sw=[30, 40.9375, -124.6875, -106.8125]) # S is 25.1875
    
def plot2map(shape,m,mapbox,color="blue",thick=1.0,label=None):
    """plot a set of polygons on a map

    INPUT:
        shape comes from shapefile.py output
        m is a basemap
        mapbox is the bounding box for the map (South,North,West,East)
        color is the line color to draw
        thick is the line width to draw
        label is the label to associate with these polygons for the legend
    OUTPUT:
        NONE
    SIDE EFFECT:
        shapes drawn on a map
    """
    plotlabel=True
    for poly in shape:
        # pull out the important parts of the shape
        # bounding box
        bbox=poly.shape.bbox
        # list of polygons within the list of points for this shape
        parts=poly.shape.parts
        # all points in this shape
        points=poly.shape.points
        
        # if the shape is not in the maps bounding box, don't draw anything
        if (bbox[2]<mapbox[2]) or (bbox[0]>mapbox[3]) or (bbox[3]<mapbox[0]) or (bbox[1]>mapbox[1]):
            pass
        else:
            # loop through all polygons within this shape
            for i in range(len(parts)):
                # set up the array of points to plot for this polygon
                if i==len(parts)-1:
                    # if this is the last part there is no next part to use as an 
                    # index so index to the end of the array
                    xypoints=np.array(points[parts[i]:])
                else:
                    # for all other points
                    xypoints=np.array(points[parts[i]:parts[i+1]])
                    
                # plot the polygon
                # first time through we need to add a label for the legend (could also be if i==0)
                if plotlabel:
                    print(label)
                    m.plot(xypoints[:,0],xypoints[:,1],color=color,linewidth=thick,label=label.replace("0","-"))
                    plotlabel=False
                else:
                    m.plot(xypoints[:,0],xypoints[:,1],color=color,linewidth=thick)
            # note, to draw shaded (filled) polygons you will need something like:
            # possibly after converting xypoints from latlon to window plot coordinates (?) with m(x,y)
            # poly = plt.Polygon(zip(xypoints[:,0],xypoints[:,1]),facecolor=colorval,edgecolor='none')
            # plt.gca().add_patch(poly)
def plot_box(map,bounds):
    x=[bounds[2],bounds[2],bounds[3],bounds[3],bounds[2]]
    y=[bounds[0],bounds[1],bounds[1],bounds[0],bounds[0]]
    map.plot(x,y,linewidth=1.3,color="black")

def shade_box(map,bounds):
    from matplotlib.patches import Rectangle
    x0, y0 = bounds[2], bounds[0]
    dx=bounds[3]-bounds[2]
    dy=bounds[1]-bounds[0]
    
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((x0,y0), dx, dy, facecolor="lightgrey",alpha=0.5))

def main():
    files=glob.glob(huc_location+"HUC*/*.shp")
    # all shape files to plot
    files.sort()
    # just use the first three
    files=files[:3]
    
    conus=True
    plot_subdomains=conus
    # Colors are speced for 5 huc levels (not just three)
    colors=["red","blue","green","magenta","cyan"]
    thicknesses=[2.5,2,1.5,1,0.5]
    # headwaters geographic bounds [South,North,West,East]
    geo=[35,43,-113,-101]
    if conus:
        # conus geographic bounds
        geo=[25,52,-125,-65]
        # and line widths
        thicknesses=[0.5,0.3,0.2]
        # and figure size
        plt.figure(figsize=(10,5))
    # set up a basemap projection to use (cylindrical) with those geographic bounds
    m = Basemap(projection='cyl',llcrnrlat=geo[0],urcrnrlat=geo[1],\
                llcrnrlon=geo[2],urcrnrlon=geo[3],resolution="i")
                
    # loop through shapefiles ploting on the map
    for i in range(len(files)):
        # position in files array isn't just in order for various (arbitrary) reasons. 
        curpos=len(files)-i-1
        filename=files[curpos]
        print(filename)
        # load the shapefile
        shp=shapefile.Reader(filename)
        # read the polygons from the shapefile
        polygons=shp.shapeRecords()
        # plot polygons on the map
        plot2map(polygons,m,geo,color=colors[curpos],
                 thick=thicknesses[curpos],
                 label=filename.split(huc_location)[-1].split("/")[0])
        
    # finally, draw ancillary stuff (grid lines, state, country, continent borders)
    if conus:
        fs=12
        m.drawparallels(np.arange(25,55,5.),labels=[1,0,0,0],dashes=[1,4],fontsize=fs)
        m.drawmeridians(np.arange(-120,-65,10.),labels=[0,0,0,1],dashes=[1,4],fontsize=fs)
        m.drawstates(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawcoastlines(linewidth=0.5)
    else:
        # for the subdomain, don't bother with country and continent borders
        m.drawstates(linewidth=1.5,color='k')
        m.drawparallels(np.arange(36,43,2.),labels=[1,0,0,0],dashes=[1,4])
        m.drawmeridians(np.arange(-112,-103,4.),labels=[0,0,0,1],dashes=[1,4])
    
    if plot_subdomains:
        for k in geosubsets.keys():
            print(k)
            shade_box(m,geosubsets[k])
        for k in geosubsets.keys():
            plot_box(m,geosubsets[k])
            
    # draw a legend
    plt.legend(loc=4)
    
    # save the file to a png
    if conus:
        region="conus"
    else:
        region="subdomain"
    plt.savefig(region+"_hucmap.png",dpi=200)


if __name__ == '__main__':
    main()