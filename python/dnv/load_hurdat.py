#!/usr/bin/env python
import datetime
import numpy as np

from bunch import Bunch

magnitudes={"WV":0,"DB":0,"LO":0,"SD":1,"TD":1,"SS":2,"TS":2,"HU":3,"EX":-1}
direction={"W":-1,"E":1,"N":1,"S":-1}

def parse_date(data):
    """Parse a datetime structure out of the input data"""
    year    =int(data[0][:4])
    month   =int(data[0][4:6])
    day     =int(data[0][6:])
    hour    =int(data[1].strip()[:2])
    minute  =int(data[1].strip()[2:])
    
    return datetime.datetime(year, month, day, hour, minute)


def add_row(data,storm):
    """Parse the current row of data and update the storm data"""
    
    strength=magnitudes[data[3].strip()]
    # find latitude and longitude, "direction" makes them negative for Western and Southern hemispheres
    lat = direction[data[4][-1]] * float(data[4][:-1])
    lon = direction[data[5][-1]] * float(data[5][:-1])
    
    # get a datetime data structure
    date = parse_date(data)
    
    # convert to km
    # 1.852 km / nautical mile
    r34 = max(0, 1.852 * np.mean(np.array(data[ 8:12], dtype=float)) )
    r50 = max(0, 1.852 * np.mean(np.array(data[12:16], dtype=float)) )
    r64 = max(0, 1.852 * np.mean(np.array(data[16:20], dtype=float)) )
    
    pressure = float(data[7])       # hPa
    wind = float(data[6]) * 0.514   # convert knots to m/s
    
    current_data=np.array(data[6:-1],dtype=int)
    
    # store the output back in the storm data structure
    storm.dates.append( date )
    storm.data.append( current_data )
    storm.lats.append( lat )
    storm.lons.append( lon )
    storm.p.append( pressure )
    storm.winds.append( wind )
    storm.r34.append( r34 )
    storm.r50.append( r50 )
    storm.r64.append( r64 )
    storm.type.append( strength )
    storm.strength = max( storm.strength, strength )
    
def new_row(data):
    return Bunch(name  = data[1].strip(),
                 strength = -2, 
                 dates = [], 
                 data  = [],
                 lats  = [],
                 lons  = [],
                 p     = [],
                 winds = [],
                 type  = [],
                 r34   = [],
                 r50   = [],
                 r64   = [])

def load(filename, year=None):
    
    outputdata=[]
    
    with open(filename,"r") as f:
        for l in f:
            data=l.split(",")
            if len(data)<5:
                if len(outputdata)>0:
                    outputdata[-1].data = np.array(outputdata[-1].data)
                outputdata.append(new_row(data))
            else:
                if len(outputdata)>0:
                    add_row(data,outputdata[-1])
                else:
                    print("ERROR: no valid rows found yet")
                    print(l)
    
    return outputdata

startdates = {
    2002:datetime.datetime(2002,10,1,0,0,0),
    2003:datetime.datetime(2003, 7,1,0,0,0),
    2004:datetime.datetime(2004, 8,1,0,0,0),
    2005:datetime.datetime(2005, 7,1,0,0,0),
    2007:datetime.datetime(2007, 9,1,0,0,0),
    2008:datetime.datetime(2008, 7,1,0,0,0)}

enddates = {
    2002:datetime.datetime(2002,10,31,23,00,00),
    2003:datetime.datetime(2003, 9,30,23,00,00),
    2004:datetime.datetime(2004, 9,30,23,00,00),
    2005:datetime.datetime(2005,10,31,23,00,00),
    2007:datetime.datetime(2007, 9,30,23,00,00),
    2008:datetime.datetime(2008, 9,30,23,00,00)}



if __name__ == '__main__':
    # 
    # sample script to plot output
    # 
    import matplotlib.pyplot as plt
    from mpl_toolkits import basemap
    
    filename="hurdat2-2000-2014-060415.txt"
    outputfile="trackmap_{}.png"
    for testyear in enddates.keys():
        plt.clf()
        print(testyear)
        startdate = startdates[testyear]
        enddate   = enddates[testyear]
        
        print("Loading data")
        data=load(filename)
        
        print("Building Map")
        m=basemap.Basemap(llcrnrlon=-100, llcrnrlat=10, urcrnrlon=-65, urcrnrlat=45, resolution="i")
        m.bluemarble(alpha=0.5)
        m.drawcoastlines(color="grey")
        m.drawmeridians(np.arange(0,-100,-5), dashes=[5,5], labels=[False, False, False,  True])
        m.drawparallels(np.arange(0,  90, 5), dashes=[5,5], labels=[ True, False, False, False])
        
        print("Plotting tracks")
        for i in data:
            if i.strength>2 and (i.dates[0]>=startdate and i.dates[-1]<=enddate):
                m.plot(i.lons,i.lats,"--",color="lightgrey")
                m.scatter(i.lons,i.lats,c=i.winds,cmap=plt.cm.jet, s=i.r64*10)
        
        plt.savefig(outputfile.format(testyear))
