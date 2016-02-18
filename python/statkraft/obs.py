from __future__ import absolute_import, print_function, division

from mpl_toolkits.basemap import pyproj
from mpl_toolkits.basemap import Basemap
import numpy as np

import mygis
from bunch import Bunch

#  assumed data file format
# double X(X) ;
#         X:units = "meters" ;
#         X:standard_name = "projection_x_coordinate" ;
#         X:long_name = "x coordinate of projection" ;
# double Y(Y) ;
#         Y:units = "meters" ;
#         Y:standard_name = "projection_y_coordinate" ;
#         Y:long_name = "y coordinate of projection" ;
# double time(time) ;
#         time:units = "days since 1900-01-01 00:00:00" ;
#         time:axis = "T" ;
#         time:calendar = "standard" ;
#         time:long_name = "time" ;
# float precipitation_amount(time, Y, X) ;
#         precipitation_amount:units = "millimeter" ;
#         precipitation_amount:missing_value = -999.99f ;
#         precipitation_amount:grid_mapping = "UTM_Zone_33" ;
#         precipitation_amount:long_name = "daily precipitation sum" ;
#         precipitation_amount:_FillValue = -999.99f ;
#         precipitation_amount:version = "1.0" ;
#         precipitation_amount:prod_date = "2015-05-15" ;
# double dummy(dummy) ;
#         dummy:units = "" ;
# double UTM_Zone_33(dummy) ;
#         UTM_Zone_33:units = "" ;
#         UTM_Zone_33:missing_value = -1. ;
#         UTM_Zone_33:grid_mapping_name = "transverse_mercator" ;
#         UTM_Zone_33:utm_zone_number = 33 ;
#         UTM_Zone_33:inverse_flattening = 298.257222101 ;
#         UTM_Zone_33:semi_major_axis = 6378137. ;
#         UTM_Zone_33:proj4 = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" ;
#         UTM_Zone_33:_CoordinateTransformType = "Projection" ;
#         UTM_Zone_33:_CoordinateAxisType = "GeoX GeoY" ;


def geo(filename, subset=False):
    X_utm=mygis.read_nc(filename,"X").data
    Y_utm=mygis.read_nc(filename,"Y").data

    X_utm,Y_utm=np.meshgrid(X_utm, Y_utm)

    projection = pyproj.Proj(proj='utm',zone=33,ellps='WGS84')
    lons, lats = projection(X_utm,Y_utm,inverse=True)

    if subset:
        lons=lons[1549:800:-1,:600]
        lats=lats[1549:800:-1,:600]

    return Bunch(lon=lons,lat=lats)
