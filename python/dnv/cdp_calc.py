#!/usr/bin/env python
import glob, datetime
import numpy as np
import load_data

input_dir = "/glade/scratch/kyoko/DNV/tracker_output/{}/TC_NATL_NEW/"
outputfile= "/glade/scratch/gutmann/dnv/cdp/{}.txt"

# CDP = 4 * ( (vm/65)**3 + 5*(Rh/50)) / vt
# for vm>65 force
#   vt>=5 
#   CDP<=10

# input file variables = 
#  YMD H lat lon slp w850 w1000 w10 vor T-T shear depth  B   VTL  VTU  rmax r64
#   0  1  2   3   4    5     6   7   8   9    10    11   12   13   14   15   16
#          [0-360]             [m/s]                                        [km]

# output file variables = 
# YMDH lat    lon     CDP anythingelse
#         [-180-180]  

def read_date(data):
    """docstring for read_date"""
    year = int(str(data[0])[:4])
    month= int(str(data[0])[4:6])
    day  = int(str(data[0])[6:8])
    hour = data[1]
    
    return datetime.datetime(year, month, day, hour)

def main(ens_member):
    """docstring for main"""
    files=glob.glob(input_dir+"track-*.txt")
    output_file=outputfile.format(ens_member)
    
    with open(output_file,"w") as output:
        for f in files:
            data=load_data.cols(f)
            npoints=data.shape[0]
            for point in range(npoints):
                if (data[point,-1] != -999):
                    startpt=max(point-1,0)
                    endpt=min(point+1,npoints-1)
                    
                    startdate = read_date(data[startpt])
                    enddate   = read_date(data[endpt])
                    dt = (enddate-startdate).total_seconds() / 3600.0 # convert seconds to hours
                    
                    distance = calc_distance(data[startpt],data[endpt])
                    
                    translation_speed = distance / (dt * 24) 
                    vt=translation_speed
                    
                    Rh=data[point,r_column]
                    vm = data[point,vmax_col]
                    cdp = 4.0 * ( (vm/65.0)**3.0 + 5.0*(Rh/50.0)) / vt
                    
                    output.write("{}  {}  {}  {}".format(
                                   str(data[point,datecol])+str(data[point,timecol]),
                                   data[point,latcol], data[point,loncol], cdp ))

if __name__ == '__main__':
    member=sys.argv[1]
    global input_dir
    input_dir.format(member)
    
    main(member)