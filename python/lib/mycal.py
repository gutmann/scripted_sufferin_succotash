# handful of routines to help managing calendars, including gregorian, noleap, 360-day etc. 
import numpy as np
import mygis
from netcdftime import utime,datetime
from bunch import Bunch



standard_calendar_names=dict(standard  = "gregorian",
                             gregorian = "gregorian",
                             noleap    = "noleap")

# names beginning with a number can't be added to a dict in the initialization (that I know of)
numerical_names=["365_day","365day","360_day","360day"]
standard_names =[ "noleap","noleap", "360day","360day"]
for i,n in zip(numerical_names,standard_names):
    standard_calendar_names[i]=n

def convert_character_dates(data):
    """Converts character date arrays (like WRF outputs) into days since 1950-01-01"""
    
    units="days since 1950-01-01 00:00:00"
    try:
        calendar=data.calendar
    except:
        calendar="gregorian"
    
    cdftime=utime(units,calendar)
    dates=np.zeros(data.shape[0])
    for i in range(data.shape[0]):
        date_string="".join(data[i,:])
        date_string=date_string.replace(":"," ").replace("/"," ").replace("-"," ").replace("_"," ")
        date_list=date_string.split()
        dateparts=[int(date_part) for date_part in date_list]
        
        dates[i] = cdftime.date2num(datetime(*dateparts))
        
    return Bunch(data=dates,calendar=calendar,units=units)

def load_time_var(filename):
    """Try to load the time netcdf variable testing a variety of names"""
    possible_time_names=["time","times","Time","Times"]
    for time_var in possible_time_names:
        try:
            time_data=mygis.read_nc(filename,time_var,returnNCvar=True)
            break
        except KeyError:
            time_data=None
    
    if time_data==None:
        raise KeyError("Time variable not found in :"+filename)
    
    if time_data.data.dtype==np.dtype('|S1'):
        output_data=convert_character_dates(time_data.data)
    else:
        output_data=Bunch(data=time_data.data[:],
                          calendar=time_data.data.calendar,
                          units=time_data.data.units)
        
    time_data.ncfile.close()
    return output_data
    
# def load_calendar(time_data):
#     try:
#         calendar=time_data.calendar
#     except:
#         calendar="gregorian"
#     return standard_calendar_names[calendar]
#
# def load_base_date(time_data):
#     """read the units attribute to find the reference/base date"""
#     units_att=time_data.units
#     if re.match("days since.*",units_att):
#         date_string=" ".join(units_att.split()[2:])
#     else:
#         date_string=units_att
#
#     date_string=date_string.replace(":"," ").replace("/"," ").replace("-"," ")
#     date_list=date_string.split()
#     dateparts=[int(date_part) for date_part in date_list)]
#
#     # passes as many parts as we have to datetime e.g. year month day or y,m,d,hour minute
#     return datetime.datetime(*dateparts)
#
# def create_dates(data,calendar,base_date):
#     """docstring for create_dates"""
#     from netcdftime import utime
#     ... just use utime from the beginning

def read_times(filename):
    """Read times from a netcdf file and convert to datetime objects"""
    
    time_data=load_time_var(filename)
    
    cdftime = utime(time_data.units,time_data.calendar)
    
    dates   = [cdftime.num2date(time_data.data[i]) for i in range(time_data.data.shape[0])]
    
    return dates
    
