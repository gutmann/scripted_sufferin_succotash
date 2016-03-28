#!/usr/bin/env python
import glob, datetime, sys
import numpy as np
import load_data

from mpl_toolkits.basemap import Basemap

global input_dir
basin="atlc"
# basin="epac"
if basin=="atlc":
    basin_subdir="TC_NATL_NEW2/"
else:
    basin_subdir="TC_EPAC_NEW2/"
# input_dir = "/glade/scratch/kyoko/DNV/tracker_output/{}/"+basin_subdir
input_dir = "/glade/scratch/kyoko/DNV/tracker_output_pgw/{}/"+basin_subdir
# input_dir = "/scratch/WEEKLY/kyoko/tracker_output_curr/{}/TC_NATL_NEW/"
# input_dir = "/scratch/WEEKLY/kyoko/tracker_output_pgw/{}/TC_NATL_NEW/"
outputfile= "cdp_ebr64_"+basin+"_{}.txt"
km_to_naut_miles = 0.539957

datecol=0
timecol=1
latcol=2
loncol=3
vmax_col=7
r_column=16


# Haversine method to calculate distance between two lat/lon points
# given by lons, lats 
def get_dist(lon1,lon2,lat1,lat2): 
    # great circle distance. 
    arg = np.sin(lat1)*np.sin(lat2)+np.cos(lat1)*np.cos(lat2)*np.cos(lon1-lon2) 
    arg = np.where(np.fabs(arg) < 1., arg, 0.999999) 
    return np.arccos(arg) 


# input file variables = 
#  YMD H lat lon slp w850 w1000 w10 vor T-T shear depth  B   VTL  VTU  rmax r64
#   0  1  2   3   4    5     6   7   8   9    10    11   12   13   14   15   16
#          [0-360]             [m/s]                                        [km]

# output file variables = 
# YMDH lat    lon     CDP anythingelse
#         [-180-180]  

def load_projection():
    """load a basemap instance for the current projection"""
    return Basemap(width=12000000,height=9000000,projection='lcc',
                resolution='c',lat_1=30.,lat_2=60,lat_0=25,lon_0=-91.)
    # return Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
    #             llcrnrlon=-180,urcrnrlon=180,lat_ts=30,resolution='c')

def read_date(data):
    """docstring for read_date"""
    year = int(str(data[datecol])[:4])
    month= int(str(data[datecol])[4:6])
    day  = int(str(data[datecol])[6:8])
    hour = int(data[timecol])
    
    return datetime.datetime(year, month, day, hour)
    
def calc_distance(startdata, enddata, projection=None):
    """docstring for calc_distance"""
    lat0,lon0 = startdata[latcol], 360 - startdata[loncol]
    lat1,lon1 =   enddata[latcol], 360 -   enddata[loncol]
    
    if (projection==None):
        projection=load_projection()
        
    x0,y0 = projection(lon0,lat0)
    x1,y1 = projection(lon1,lat1)
    
    return np.sqrt((x0-x1)**2 + (y0-y1)**2) / 1000.0 # convert m to km
    

def main(ens_member):
    """docstring for main"""
    files=glob.glob(input_dir+"track-*.txt")
    output_file=outputfile.format(ens_member)
    print(output_file)
    projection=load_projection()
    with open(output_file,"w") as output:
        for f in files:
            data=load_data.cols(f)
            npoints=data.shape[0]
            for point in range(npoints):
                if (data[point,-1] > -990):
                    startpt=max(point-1,0)
                    endpt=min(point+1,npoints-1)
                    
                    startdate = read_date(data[startpt])
                    enddate   = read_date(data[endpt])
                    # time between points in hours
                    dt = (enddate-startdate).total_seconds() / 3600.0 # convert seconds to hours
                    
                    # distance in nautical miles
                    distance = calc_distance(data[startpt],data[endpt], projection=projection) * km_to_naut_miles
                    
                    # translation speed in knots (nautical miles per hour)
                    translation_speed = distance / dt
                    vt=translation_speed
                    
                    # R64 for the current point
                    Rh=data[point,r_column]
                    vm = data[point,vmax_col]
                    
                    # CDP = 4 * ( (vm/65)**3 + 5*(Rh/50)) / vt
                    # for vm>65 force
                    #   vt>=5 
                    #   CDP<=10
                    vt=max(0.01,vt)
                    if (vm>65):vt=max(5,vt)
                    cdp = 4.0 * ( (vm/65.0)**3.0 + 5.0*(Rh/50.0)) / vt
                    cdp = min(cdp,10)
                    
                    output.write("{}  {}  {}  {}\n".format(
                                   str(int(data[point,datecol]))+str(100+int(data[point,timecol]))[1:],
                                   data[point,latcol], data[point,loncol], cdp ))

if __name__ == '__main__':
    member=sys.argv[1]
    input_dir = input_dir.format(member)
    
    main(member)