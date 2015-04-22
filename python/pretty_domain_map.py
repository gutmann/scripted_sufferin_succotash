#!/usr/bin/env python

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from myshade import shade

import sys

# Command line options.
if len(sys.argv)>2:
    center_lat=float(sys.argv[1])
    center_lon=float(sys.argv[2])
else:
    center_lat=40.
    center_lon=-105.

if len(sys.argv)>4:
    width=float(sys.argv[3])
    height=float(sys.argv[4])
else:
    width=1000.0 #km
    height=1000.0 #km

if len(sys.argv)>5:
    titlename=sys.argv[5]
else:
    titlename=""    

if len(sys.argv)>6:
    rescale_shading=(sys.argv[6]=="rescale_shading")
else:
    rescale_shading=False

if len(sys.argv)>7:
    large=(sys.argv[7]=="large_output")
else:
    large=False


# rescale_shading=True
print(titlename)

topo_dataset="/glade/u/home/gutmann/work/topo_1km_global.nc"
if (rescale_shading==False):
    rescale_colors=True
else:
    rescale_colors=False

print(large)
if large:
    print("Making a larger figure")
    plt.figure(figsize=(10,8))
else:
    plt.figure(figsize=(5,5))
title_font_size=10
axis_font_size=10
width_in_meters=width*1000.0
height_in_meters=height*1000.0

lat1=center_lat+height/333.0
lat2=center_lat-height/333.0

# distance between grid lines
def calc_grid_lines(length_scale):
    """docstring for calc_grid_lines"""
    if length_scale<=50:   return 0.1
    if length_scale<=100:  return 0.25
    if length_scale<=200:  return 0.5
    if length_scale<=400:  return 1.0
    if length_scale<=1000: return 2.0
    return 5.0


dlat=calc_grid_lines(height)
dlon=calc_grid_lines(width)
# dlat=0.25
# dlon=0.25

# setup Lambert Conformal basemap.
print("Setting up map projection")
m = Basemap(width=width_in_meters,height=height_in_meters,projection='lcc',
            resolution='h',lat_1=lat1,lat_2=lat2,lat_0=center_lat,lon_0=center_lon)

data=Dataset(topo_dataset)
print("Reading topo coordinates")
lons=data.variables["X"][5000:15000]
lats=data.variables["Y"][10000:0:-1]
# if the area to be displayed doesn't fit in the subdomain, we have to read in more data
if (   ((center_lat + height/222) > lats.max())
    or ((center_lat - height/222) < lats.min())
    or ((center_lon + width/222)  > lons.max())
    or ((center_lon - width/222)  < lons.min())):
    lons=data.variables["X"][:]
    lats=data.variables["Y"][::-1]
    print("Reading FULL topo data")
    topoin=data.variables["topo"][::-1,:]
else:
    print("Reading topo data")
    topoin=data.variables["topo"][10000:0:-1,5000:15000]

print("Setting up plot data")
nx = int((m.xmax-m.xmin)/1000.)+1; ny = int((m.ymax-m.ymin)/1000.)+1

# transform to nx x ny regularly spaced 1km native projection grid
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# if type(topodat)==np.ma.core.MaskedArray:topodat.mask=False
zc=topodat.copy()
topodat[topodat<0]=0

# create the pretty color shaded map
if rescale_colors:
    colormax=min(5000,topodat.max())
    midpoint=(colormax+topodat.min())/2
    colorrange=(colormax-topodat.min()) * 0.90
    colormin=midpoint-colorrange
    colormax=midpoint+colorrange/2
else:
    toposcalar=1.0
    if topodat.max()<2500:toposcalar=0.6
    colormin=-1800 * toposcalar
    colormax=3500  * toposcalar

vertical_exageration=1.0
if rescale_shading:
    topomax=topodat.max()
    toporange=topomax-topodat.min()
    # topomin=topomax-(toporange/0.64)
    vertical_exageration=min(max(1000.0/toporange,1),5)
    print("Vertical Exageration = "+str(vertical_exageration))


zc[zc<=0]=colormin
img=shade(topodat*vertical_exageration,colordata=zc,clim=(colormin,colormax),dx=185,cmap=plt.cm.terrain)
# make oceans transparent
img[:,:,3][topodat<=0]=0

print("Drawing map")
# draw a boundary around the map, fill the background.
# this background will end up being the ocean color, since
# the continents will be drawn on top.
m.drawmapboundary(fill_color='lightblue')

# fill continents, set lake color same as ocean color.
m.fillcontinents(lake_color="blue",color=(0,0,0,0));

map_img=m.imshow(img,origin="lower")

if rescale_colors:
    map_img.set_cmap("terrain")
    map_img.set_clim(colormin,colormax)
    m.colorbar()

# draw parallels and meridians.
# label parallels on right and left
# label meridians on bottom
parallels = np.arange(-88.,88,dlat)
# labels = [left,right,top,bottom]
if rescale_colors:
    m.drawparallels(parallels,labels=[True,False,False,False],fontsize=axis_font_size)
else:
    m.drawparallels(parallels,labels=[True,True,False,False],fontsize=axis_font_size)
meridians = np.arange(0.,358.,dlon)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=axis_font_size)

if large:
    m.drawstates(linewidth=1.5)
    m.drawcountries(linewidth=2)
    # m.drawcoastlines(linewidth=2)
    m.drawrivers(color="Blue",linewidth=0.5)
else:
    m.drawstates(linewidth=1.5)
    m.drawcountries(linewidth=2)
    m.drawrivers(color="Blue",linewidth=1)
    
plt.title(titlename,fontsize=title_font_size)

xpt,ypt = m(center_lon,center_lat)
# if not large:
    # m.plot(xpt,ypt,'ro')  # plot a red dot at the center point
# put some text next to the dot, offset a little bit
# (the offset is in map projection coordinates)
# plt.text(xpt+100000,ypt+100000,'Boulder (%5.1fW,%3.1fN)' % (lonpt,latpt))

try:
    ID=titlename.split("(")[1]
    ID=ID.replace(")","").replace(" ","")
except:
    ID=titlename
if not ID:
    ID="test"

print(ID+".dem.png")
plt.savefig(ID+".dem.png",dpi=200)


# plot blue dot on e.g. Boulder, Colorado and label it as such.
# lon, lat = -105.278, 40.019 # Location of Boulder
# convert to map projection coords.
# Note that lon,lat can be scalars, lists or numpy arrays.
# xpt,ypt = m(lon,lat)
# m.plot(xpt,ypt,'bo')  # plot a blue dot there
# put some text next to the dot, offset a little bit
# (the offset is in map projection coordinates)
# plt.text(xpt+100000,ypt+100000,'Boulder (%5.1fW,%3.1fN)' % (lonpt,latpt))
# plt.show()