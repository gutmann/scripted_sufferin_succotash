#!/usr/bin/env python

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from myshade import shade
import custom_cmap

import sys

def usage():
    """docstring for usage"""
    print("""
pretty_domain_maps.py center_lat center_lon width height name rescale_shading large_output
    
    center_lat = latitude coordinate of map center  [default=40]
    center_lon = longitude coordinate of map center [default=-105]
    width      = domain width E-W (km)              [default=1000]
    height     = domain height N-S (km)             [default=1000]
    name       = map title (and filename)           [default='']
    rescale_shading     Boolean flag if present shading is rescaled for local topography
    large_output        Boolean flag if present a larger output image is created
    
    """)

# Command line options.
if len(sys.argv)>1:
    if sys.argv[1]=="-h" or sys.argv[1]=="--help":
        usage()
        sys.exit()
else:
    usage()
    sys.exit()

if len(sys.argv)>2:
    center_lat=float(sys.argv[1])
    center_lon=float(sys.argv[2])
else:
    center_lat=39.
    center_lon=-107.

if len(sys.argv)>4:
    width=float(sys.argv[3])
    height=float(sys.argv[4])
else:
    width=1268.0 #km
    height=1052.0 #km

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

print("Creating a 'large' map="+str(large))
rescale_colors=False
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
    if length_scale<=1200: return 2.0
    return 4.0


dlat=calc_grid_lines(height)
dlon=calc_grid_lines(width)
# dlat=0.25
# dlon=0.25
dx=2000.0

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
nx = int((m.xmax-m.xmin)/dx)+1; ny = int((m.ymax-m.ymin)/dx)+1

# transform to nx x ny regularly spaced dx O(1-4km) native projection grid
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
    colormin=500 # 0 # -1800 * toposcalar
    colormax=3800  * toposcalar

vertical_exageration=1.0
if rescale_shading:
    topomax=topodat.max()
    toporange=topomax-topodat.min()
    # topomin=topomax-(toporange/0.64)
    vertical_exageration=min(max(1000.0/toporange,1),5)
    print("Vertical Exageration = "+str(vertical_exageration))


zc[zc<=colormin+20]=colormin+20
img=shade(topodat*vertical_exageration,colordata=zc,clim=(colormin,colormax),dx=185,cmap=custom_cmap.terrain())
# make oceans transparent
img[:,:,3][topodat<=0]=0

print("Drawing map")
# draw a boundary around the map, fill the background.
# this background will end up being the ocean color, since
# the continents will be drawn on top.
m.drawmapboundary(fill_color='lightblue')

# fill continents, set lake color
m.fillcontinents(lake_color="blue",color=(0,0,0,0));

map_img=m.imshow(img,origin="lower")

# print("rescale_colors="+str(rescale_colors))
if rescale_colors:
    map_img.set_cmap(custom_cmap.terrain()) #"terrain")
    map_img.set_clim(colormin,colormax)
    m.colorbar()
# else:
#     map_img.set_cmap(custom_cmap.terrain()) #"terrain")
#     map_img.set_clim(colormin,colormax)
#     m.colorbar()

# draw parallels and meridians.
# label parallels on right and left
# label meridians on bottom
parallels = np.arange(-88.,88,dlat)
# labels = [left,right,top,bottom]
if rescale_colors:
    m.drawparallels(parallels,linewidth=0.5, labels=[True,False,False,False],fontsize=axis_font_size)
else:
    m.drawparallels(parallels,linewidth=0.5, labels=[True,False,False,False],fontsize=axis_font_size)
meridians = np.arange(0.,358.,dlon)
m.drawmeridians(meridians,linewidth=0.5, labels=[False,False,False,True],fontsize=axis_font_size)

if large:
    m.drawstates(linewidth=1.5)
    m.drawcountries(linewidth=2)
    # m.drawcoastlines(linewidth=2)
    m.drawrivers(color="blue",linewidth=0.5)
else:
    m.drawstates(linewidth=1)
    m.drawcountries(linewidth=2)
    m.drawrivers(color="blue",linewidth=0.2)
    
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
plt.savefig(ID+".dem.png",dpi=300)


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