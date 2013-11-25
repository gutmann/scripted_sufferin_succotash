#!/usr/bin/env python
import glob,os,shutil
import numpy as np

import swim_io as io
import units
from bunch import Bunch

wrf_directory="/glade/p/ral/RHAP/asd000/HW2010.2/"

lat=40.000
lon=-105.1

template_file="template.dat"
output_file="met_data_"+str(lat)+"_"+str(lon)".dat"
wrf_file_search=wrf_directory+"ctrl/wrfout/wrfout*"
wrf_geo_grid_file=wrf_directory+"geo_em.d01-4km.nc"

searchstrings=["__STARTDATE__","__ENDDATE__","__LATITUDE__","__LONGITUDE__",
               "__TSOIL1__","__TSOIL2__","__TSOIL3__","__TSOIL4__",
               "__MSOIL1__","__MSOIL2__","__MSOIL3__","__MSOIL4__",
               "__LSOIL1__","__LSOIL2__","__LSOIL3__","__LSOIL4__",
               "__SKINTEMP__","__DEEPSOILTEMP__",
               "__VEGTYPE__","__SOILTYPE__","__SLOPETYPE__",
               "__VEGFRAC__"]

    # ---------------------------------------------------
    # VARIABLES TO READ FROM geo_grid_file
    # ---------------------------------------------------
    # float XLAT_M(Time, south_north, west_east) ;
    # float XLONG_M(Time, south_north, west_east) ;
    # 
    # float GREENFRAC(Time, month, south_north, west_east) ;
    # float SLOPECAT(Time, south_north, west_east) ;
    
    # ---------------------------------------------------
    # VARIABLES TO READ FROM wrf output files
    # ---------------------------------------------------
    # Position / Time data
    # ---------------------------------------------------
    # XLAT,XLONG  float dims(Time, south_north, west_east) ;
    # Times char dims(Time,19) e.g. "2008-12-28_00:00:00"
    
    # ---------------------------------------------------
    # Met Forcing data
    # ---------------------------------------------------
    # float dims(Time, south_north, west_east)
    # Q2:description = "QV at 2 M" ; kg kg-1
    # T2:description = "TEMP at 2 M" ; K
    # PSFC:description = "SFC PRESSURE" ; Pa
    # U10:description = "U at 10 M" ; m/s
    # V10:description = "V at 10 M" ; m/s
    # RAINNC:description = "ACCUMULATED TOTAL GRID SCALE PRECIPITATION" ; mm
    # SWDOWN:description = "DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE" ; W/m2
    # GLW:description = "DOWNWARD LONG WAVE FLUX AT GROUND SURFACE" ; W/m2
    
    # ---------------------------------------------------
    # Essentially State data
    # ---------------------------------------------------
    # TSK:description = "SURFACE SKIN TEMPERATURE" ; K
    # VEGFRA:description = "VEGETATION FRACTION" ; [m2/m2]
    
    # float dims(Time, soil_layers_stag, south_north, west_east)
    # TSLB:description = "SOIL TEMPERATURE" ; K
    # SMOIS:description = "SOIL MOISTURE" ; m3/m3
    # SH2O:description = "SOIL LIQUID WATER" ; m3/m3
    
    # int dims(Time, south_north, west_east)
    # IVGTYP:description = "DOMINANT VEGETATION CATEGORY" ; []
    # ISLTYP:description = "DOMINANT SOIL CATEGORY" ; []
    
    
    
#   DATE           sqrt(U10^2+V10^2)        0                T2          f(Q2,PSFC,T2)      PSFC/100           SWDOWN          GLW             f(RAINNC)
#  -----------------------------------------------------------------------------------------------------------------------------------------------
#   UTC date/time        windspeed       wind dir         temperature      humidity        pressure           shortwave      longwave       precipitation
#  yyyy mm dd hh mi       m s{-1}        degrees               K               %             hPa               W m{-2}        W m{-2}       kg m{-2} s{-1}
#  -----------------------------------------------------------------------------------------------------------------------------------------------
#  <Forcing>  This tag ("<Forcing>", not case sensitive) begins the section of forcing data.  
#  1998 01 01 06 30     5.6300001144   178.0000000000   263.9499816895    86.0999984741  1002.0000000000     0.0000000000   281.0000000000     0.0000000000
#     

def replace_in_file(filename,searchstring,outputstring):
    """replace searchstring with outputstring in the file filename
    
    Read all lines in filename, replaceing searchstring with outputstring 
    and writing to a new temporary file as you go, then copy temp back to filename
    """
    with open(filename,"ru") as f:
        with open(filename+".temporary","wu") as output:
            for l in f:
                if type(searchstring)==list:
                    for s,o in zip(searchstring,outputstring):
                        l=l.replace(s,o)
                else:
                    l=l.replace(searchstring,outputstring)
                output.write(l)
    os.rename(filename+".temporary",filename)

def date_from_filename(filename):
    """docstring for date_from_filename"""
    # filename format=wrfout_d01_2008-12-19_00:00:00.nc
    fileparts=filename.split("-")
    year=fileparts[0][-4:]
    month=fileparts[1]
    day=fileparts[2][:2]
    return year+month+day
    
def get_latlon(filename,lat,lon,latvar="XLAT",lonvar="XLONG"):
    """subset out the time dimention from lat and lon variables"""
    latdata=swim_io.read_nc(filename,latvar,returnNCvar=True)
    latgrid=latdata.data[0,...]
    latdata.ncfile.close()
    londata=swim_io.read_nc(filename,lonvar,returnNCvar=True)
    longrid=londata.data[0,...]
    londata.ncfile.close()
    
    dists=(latgrid-lat)**2 + (longrid-lon)**2
    y,x=np.unravel_index(np.argmin(dists),dists.shape)
    
    return y,x
    

def get_var_from_file_at_pos(variable,level,lat,lon,filename,latvar="XLAT",lonvar="XLONG"):
    """Read data from a file with assumptions about variable grid
    return data at matching lat,lon
    
    assumes
        dims(Time, [levels], south_north, west_east)
        
    if level==0 assume levels dim doesn't exist    
    """
    y,x=get_latlon_pos(filename,lat,lon,latvar=latvar,lonvar=lonvar)
    
    vardata=swim_io.read_nc(filename,variable,returnNCvar=True)
    if level>0:
        results=vardata.data[0,level-1,y,x]
    else:
        results=vardata.data[0,y,x]
    vardata.ncfile.close()
    return results

    
def load_parameters(geogrid_file,wrf_filesearch,lat,lon):
    """Read noahmp parameters from a wrf geogrid (netcdf) file
    
    Reads/Stores/calculates: 
                STARTDATE , ENDDATE , LATITUDE , LONGITUDE ,
                TSOIL1 , TSOIL2 , TSOIL3 , TSOIL4 ,
                MSOIL1 , MSOIL2 , MSOIL3 , MSOIL4 ,
                LSOIL1 , LSOIL2 , LSOIL3 , LSOIL4 ,
                SKINTEMP , DEEPSOILTEMP ,
                VEGTYPE , SOILTYPE , SLOPETYPE ,
                VEGFRAC
    """
    files=glob.glob(wrf_filesearch)
    files.sort()
    
    parameters=dict()
    parameters["STARTDATE"] = date_from_filename(files[0])
    parameters["ENDDATE"]   = date_from_filename(files[-1])
    parameters["LATITUDE"]  = lat
    parameters["LONGITUDE"] = lon
    for i in range(1,5):
        # TSLB:description = "SOIL TEMPERATURE" ; K
        parameters["TSOIL"+str(i)]=get_var_from_file_at_pos("TSLB",level=i,lat=lat,lon=lon,filename=files[0])
        # SMOIS:description = "SOIL MOISTURE" ; m3/m3
        parameters["MSOIL"+str(i)]=get_var_from_file_at_pos("SMOIS",level=i,lat=lat,lon=lon,filename=files[0])
        # SH2O:description = "SOIL LIQUID WATER" ; m3/m3
        parameters["LSOIL"+str(i)]=get_var_from_file_at_pos("SH2O",level=i,lat=lat,lon=lon,filename=files[0])
        
    # TSK:description = "SURFACE SKIN TEMPERATURE" ; K
    parameters["SKINTEMP"]=get_var_from_file_at_pos("TSK",level=0,lat=lat,lon=lon,filename=files[0])
    # SOILTEMP:description = "Annual mean deep soil temperature" ; K
    parameters["DEEPSOILTEMP"]=get_var_from_file_at_pos("SOILTEMP",level=0,lat=lat,lon=lon,
                                                        filename=geogrid_file,latvar="XLAT_M",lonvar="XLONG_M")
    # IVGTYP:description = "DOMINANT VEGETATION CATEGORY" ; []
    parameters["VEGTYPE"]=get_var_from_file_at_pos("IVGTYP",level=0,lat=lat,lon=lon,filename=files[0])
    # ISLTYP:description = "DOMINANT SOIL CATEGORY" ; []
    parameters["SOILTYPE"]=get_var_from_file_at_pos("ISLTYP",level=0,lat=lat,lon=lon,filename=files[0])
    # SLOPECAT:description = "Dominant category" ; []
    parameters["SLOPETYPE"]=get_var_from_file_at_pos("SLOPECAT",level=0,lat=lat,lon=lon,
                                                        filename=geogrid_file,latvar="XLAT_M",lonvar="XLONG_M")
    maxgf=0.0
    for i in range(1,13):
        gf=get_var_from_file_at_pos("GREENFRAC",level=i,lat=lat,lon=lon,
                                    filename=geogrid_file,latvar="XLAT_M",lonvar="XLONG_M")
        maxgf=max(gf,maxgf)
        
    parameters["VEGFRAC"]=maxgf
    
    return parameters
    

def get_timeseries(files,varname,y,x):
    
    # read data from first file to get the number of time steps per file 
    # (ideally this should be allowed to vary... it sort of is) but we don't test for overflow
    data=swim_io.read_nc(files][0],varname,returnNCvar=True)
    # get data at pos
    initdata=data.data[:,y,x]
    # close the netcdf file
    data.ncfile.close()
    
    # the number of time steps per file is simply the length of the current data
    timesPerFile=initdata.size
    
    # initialize the full time series of data and store the first files worth
    outputdata=np.zeros(len(files)*timesPerFile)
    outputdata[:initdata.size]=initdata
    curpos=initdata.size
    # loop over files storing the rest of the data
    for i in range(1,len(files)):
        # open file for this variable
        data=swim_io.read_nc(files][0],varname,returnNCvar=True)
        # store this time series
        outputdata[curpos:curpos+data.data.shape[0]]=data.data[:,y,x]
        curpos+=data.data.shape[0]
        # close the netcdf file
        data.ncfile.close()
        
    return outputdata[:curpos]
        
    


def read_times(files,varname):
    dates=np.concatenate(swim_io.read_files(files,varname))
    alldates=[]
    for i in range(dates.shape[0]):
        alldates.append("".join(dates[i][:16]))
        alldates[i].replace("-"," ")
        alldates[i].replace("_"," ")
        alldates[i].replace(":"," ")
    
    return alldates
    
def load_forcing_data(filesearch,lat,lon):
    """Load a time series of WRF forcing data for the lat/lon position given"""
    files=glob.glob(filesearch)
    files.sort()
    y,x=get_latlon(files[0],lat,lon)
    
    u=get_timeseries(files,"U10",y,x)
    v=get_timeseries(files,"V10",y,x)
    windspeed=np.sqrt(u**2+v**2)
    t=get_timeseries(files,"T2",y,x)
    p=get_timeseries(files,"PSFC",y,x)/100.0
    q=get_timeseries(files,"Q2",y,x)
    rh=units.mr2rh(t,p,q)*100.0
    precip=get_timeseries(files,"RAINNC",y,x)/3600.0
    precip[1:]=precip[1:]-precip[:-1]
    sw=get_timeseries(files,"SWDOWN",y,x)
    lw=get_timeseries(files,"GLW",y,x)
    
    times=read_times(file,"Times")

    return Bunch(wind=windspeed,winddir=windspeed*0, t=t,rh=rh,p=p,precip=precip,sw=sw,lw=lw,times=times)
    
def write_forcing(filename,forcing):
    infile=open(filename)
    outfile=open(filename+".temp","wu")
    for l in infile:
        outfile.write(l)
    infile.close()
    
#   UTC date/time        windspeed       wind dir         temperature      humidity        pressure           shortwave      longwave       precipitation
#  1998 01 01 06 30     5.6300001144   178.0000000000   263.9499816895    86.0999984741  1002.0000000000     0.0000000000   281.0000000000     0.0000000000
    formatstring="{0.times[_i_]  0.wind[_i_] 0.winddir[_i_] 0.t[_i_] 0.rh[_i_] 0.p[_i_] 0.sw[_i_] 0.lw[_i_] 0.precip[_i_]}"
    for i in range(forcing.wind.size):
        curformat=formatstring.replace("_i_",str(i))
        outfile.write(curformat.format(forcing))
    outfile.close
    os.rename(filename+".temp",filename)
        
    

def write_outputfile(template,outputfile,parameters,forcing):
    """docstring for write_outputfile"""
    shutil.copyfile(template,outputfile)
    for s in searchstrings:
        replace_in_file(outputfile,s,parameters[s[2:-2]])
    write_forcing(outputfile,forcing)
    

def main():
    """Write a noahmp Simple Driver input file from WRF geogrid and outputfiles"""
    parameters=load_parameters(wrf_geo_grid_file,lat,lon)
    forcing=load_forcing_data(wrf_file_search,lat,lon)
    write_outputfile(template_file, output_file, parameters,forcing)

if __name__ == '__main__':
    main()