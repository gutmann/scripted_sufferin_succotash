#!/usr/bin/env python
import datetime
import numpy as np
import mygis
import glob
from bunch import Bunch

outputfile = "summary_{}_{}_{}"
models = ["CCSM4","CNRM-CM5","GFDL-CM3","MIROC5","NorESM1-M"]
scenarios = ["historical","rcp45","rcp85"]
varnames = ["pr","tas"]

month_lengths = np.array([0,31,28,31,30,31,30,31,31,30,31,30,31])
month_starts = month_lengths.cumsum()

start_years={"historical":1980, 
             "rcp45":2020,
             "rcp85":2020}
end_years={"historical":2005,
             "rcp45":2050,
             "rcp85":2050}


def parse_noleap(times, startstring):
    init_date   = startstring.split("-")
    start_year  = int(init_date[0])
    start_month = int(init_date[1])
    start_day   = int(init_date[2][:2])
    t0 = start_year*365 + month_starts[start_month-1] + start_day
    
    dates = times + t0
    
    years = (dates / 365).astype('i')
    
    doy = dates - (years * 365)
    months = np.zeros(doy.shape,dtype=int)
    
    for i in range(1,len(month_starts)):
        is_month = np.where((month_starts[i]>doy) & (month_starts[i-1]<=doy))
        months[is_month] = i
    
    return years, months

def parse_gregorian(times, startstring):
    init_date   = startstring.split("-")
    start_year  = int(init_date[0])
    start_month = int(init_date[1])
    start_day   = int(init_date[2][:2])
    try:
        start_hour   = int(init_date[2].split(":")[0][-2:])
        start_minute = int(init_date[2].split(":")[1])
        try:
            start_second   = int(init_date[2].split(":")[2][:2])
        except:
            start_second = 0
    except:
        start_hour = 0
        start_minute = 0
        start_second = 0
        
    t0 = datetime.datetime(start_year, start_month, start_day, start_hour, start_minute, start_second)
    dates = [t0 + datetime.timedelta(t) for t in times]
    years = np.array([d.year for d in dates])
    months = np.array([d.month for d in dates])
    return years, months
    

def parse_times(times, start, calendar):
    if calendar=="noleap" or (calendar=="365day"):
        return parse_noleap(times, start)
    elif (calendar=="gregorian") or (calendar=="standard"):
        return parse_gregorian(times, start)

def read_times(filename):
    atts = mygis.read_atts(filename,"time")
    calendar = atts.calendar
    start_date = atts.units.split(" since ")[1]
    
    times = mygis.read_nc(filename,"time").data
    
    return parse_times(times, start_date, calendar)

def main(model, scenario, varname):
    files=glob.glob(model+"/"+scenario+"/day/atmos/day/r1i1p1/latest/"+varname+"/*.nc")
    files.sort()
    
    atts = mygis.read_atts(files[0],varname)
    gatts = mygis.read_atts(files[0],global_atts=True)
    geo  = mygis.read_geo(files[0])
    lat = Bunch(data=geo.lat,name="lat",dtype='f',attributes=geo.latatts,dims=('lat','lon'))
    lon = Bunch(data=geo.lon,name="lon",dtype='f',attributes=geo.lonatts,dims=('lat','lon'))
    datadims=('time','lat','lon')
    start_year = start_years[scenario]
    end_year   = end_years[scenario]
    
    print("Getting time data")
    years, months = read_times(files[0])
    print(files[0].split("/")[-1])
    print(years[0],years[-15:])
    print(months[0],months[-15:])
    
    print("Reading {} data".format(varname))
    data = mygis.read_nc(files[0],varname,returnNCvar=True)
    output = np.zeros((13,data.data.shape[1],data.data.shape[2]))
    data.ncfile.close()
    n = np.zeros(13)
    for f in files:
        print("    "+f)
        years, months = read_times(f)
        data = mygis.read_nc(f,varname,returnNCvar=True)
        for i in range(data.data.shape[0]):
            if (years[i]>start_year) & (years[i]<=end_year):
                output[months[i]-1] += data.data[i]
                n[months[i]-1]+=1
        data.ncfile.close()
    
    output[-1]=output[:-1].sum(axis=0)
    n[-1] = n[:-1].sum(axis=0)
    if varname=="pr":
        n[-1]/=365.0
        n[:-1] /= month_lengths[1:]
        n/=86400.0
    print("writing")
    results = output/n[:,np.newaxis,np.newaxis]
    mygis.write(outputfile.format(model, scenario,varname),results, dims=datadims, extravars=[lat,lon],
                varname=varname, attributes=atts, global_attributes=gatts)

if __name__ == '__main__':
    
    # main("NorESM1-M",scenarios[0],varnames[0])
    for m in models[-1:]:
        for s in scenarios:
            for v in varnames:
                print(m,s,v)
                # try:
                main(m,s,v)
                # except exception as e:
                    # print(e)
