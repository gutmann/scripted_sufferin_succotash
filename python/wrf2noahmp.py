#!/usr/bin/env python
import glob,os,shutil
import datetime

import numpy as np

import swim_io
import units
from bunch import Bunch

wrf_directory="/glade/p/ral/RHAP/asd000/HW2010.2/"

template_file="template.dat"
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
    fileparts=filename.split("/")[-1]
    fileparts=fileparts.split("-")
    year=fileparts[0][-4:]
    month=fileparts[1]
    day=fileparts[2][:2]
    return year+month+day
    
def get_latlon_pos(filename,lat,lon,latvar="XLAT",lonvar="XLONG"):
    """subset out the time dimention from lat and lon variables"""
    latdata=swim_io.read_nc(filename,latvar,returnNCvar=True)
    latgrid=latdata.data[0,...]
    latdata.ncfile.close()
    londata=swim_io.read_nc(filename,lonvar,returnNCvar=True)
    longrid=londata.data[0,...]
    londata.ncfile.close()
    
    dists=(latgrid-lat)**2 + (longrid-lon)**2
    y,x=np.unravel_index(np.argmin(dists),dists.shape)
    
    return int(y),int(x)
    

def get_var_from_file_at_pos(variable,level,lat,lon,y=None,x=None,filename="",latvar="XLAT",lonvar="XLONG"):
    """Read data from a file with assumptions about variable grid
    return data at matching lat,lon
    
    assumes
        dims(Time, [levels], south_north, west_east)
        
    if level==0 assume levels dim doesn't exist    
    """
    if y==None:
        y,x=get_latlon_pos(filename,lat,lon,latvar=latvar,lonvar=lonvar)
    
    vardata=swim_io.read_nc(filename,variable,returnNCvar=True)
    if level>0:
        outputdata=vardata.data[0,level-1,y,x]
    else:
        outputdata=vardata.data[0,y,x]
    vardata.ncfile.close()
    return outputdata

    
def load_parameters(geogrid_file,wrf_filesearch,lat,lon,y=None,x=None):
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
        parameters["TSOIL"+str(i)]=get_var_from_file_at_pos("TSLB",level=i,lat=lat,lon=lon,y=y,x=x,filename=files[0])
        # SMOIS:description = "SOIL MOISTURE" ; m3/m3
        parameters["MSOIL"+str(i)]=get_var_from_file_at_pos("SMOIS",level=i,lat=lat,lon=lon,y=y,x=x,filename=files[0])
        # SH2O:description = "SOIL LIQUID WATER" ; m3/m3
        parameters["LSOIL"+str(i)]=get_var_from_file_at_pos("SH2O",level=i,lat=lat,lon=lon,y=y,x=x,filename=files[0])
        
    # TSK:description = "SURFACE SKIN TEMPERATURE" ; K
    parameters["SKINTEMP"]=get_var_from_file_at_pos("TSK",level=0,lat=lat,lon=lon,y=y,x=x,filename=files[0])
    # SOILTEMP:description = "Annual mean deep soil temperature" ; K
    parameters["DEEPSOILTEMP"]=get_var_from_file_at_pos("SOILTEMP",level=0,lat=lat,lon=lon,y=y,x=x,
                                                        filename=geogrid_file,latvar="XLAT_M",lonvar="XLONG_M")
    # IVGTYP:description = "DOMINANT VEGETATION CATEGORY" ; []
    parameters["VEGTYPE"]=get_var_from_file_at_pos("IVGTYP",level=0,lat=lat,lon=lon,y=y,x=x,filename=files[0])
    # ISLTYP:description = "DOMINANT SOIL CATEGORY" ; []
    parameters["SOILTYPE"]=get_var_from_file_at_pos("ISLTYP",level=0,lat=lat,lon=lon,y=y,x=x,filename=files[0])
    # SLOPECAT:description = "Dominant category" ; []
    parameters["SLOPETYPE"]=get_var_from_file_at_pos("SLOPECAT",level=0,lat=lat,lon=lon,y=y,x=x,
                                                        filename=geogrid_file,latvar="XLAT_M",lonvar="XLONG_M")
    maxgf=0.0
    for i in range(1,13):
        gf=get_var_from_file_at_pos("GREENFRAC",level=i,lat=lat,lon=lon,y=y,x=x,
                                    filename=geogrid_file,latvar="XLAT_M",lonvar="XLONG_M")
        maxgf=max(gf,maxgf)
        
    parameters["VEGFRAC"]=maxgf
    
    return parameters
    

def get_timeseries(files,varname,y,x):
    
    # read data from first file to get the number of time steps per file 
    # (ideally this should be allowed to vary... it sort of is) but we don't test for overflow
    data=swim_io.read_nc(files[0],varname,returnNCvar=True)
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
        data=swim_io.read_nc(files[i],varname,returnNCvar=True)
        # store this time series
        outputdata[curpos:curpos+data.data.shape[0]]=data.data[:,y,x]
        curpos+=data.data.shape[0]
        # close the netcdf file
        data.ncfile.close()
        
    return outputdata[:curpos]

def calc_times(filesearch):
    files=glob.glob(filesearch)
    files.sort()
    t1=date_from_filename(files[0])
    d1=datetime.datetime(int(t1[:4]),int(t1[4:6]),int(t1[6:8]),0,0)
    dates=[str(d1+datetime.timedelta(i/24.0)).replace("-"," ").replace(":"," ")[:-2] for i in range(len(files)*24)]
    return dates

def read_times(filesearch,varname): #THIS IS INCREDIBLY INEFFICIENT... not sure why Nio can't read string arrays fast!
    return calc_times(filesearch)
    
    # dates=np.concatenate(swim_io.read_files(filesearch,varname))
    # alldates=[]
    # for i in range(dates.shape[0]):
    #     alldates.append("".join(dates[i][:16]))
    #     alldates[i].replace("-"," ")
    #     alldates[i].replace("_"," ")
    #     alldates[i].replace(":"," ")
    # 
    # return alldates
    
def load_forcing_data(filesearch,lat,lon,data=None):
    """Load a time series of WRF forcing data for the lat/lon position given"""
    files=glob.glob(filesearch)
    files.sort()
    y,x=get_latlon_pos(files[0],lat,lon)
    print(y,x)
    
    print("u")
    u=get_timeseries(files,"U10",y,x)
    print("v")
    v=get_timeseries(files,"V10",y,x)
    windspeed=np.sqrt(u**2+v**2)
    print("t")
    t=get_timeseries(files,"T2",y,x)
    print("p")
    p=get_timeseries(files,"PSFC",y,x)/100.0
    print("q")
    q=get_timeseries(files,"Q2",y,x)
    rh=units.mr2rh(t,p,q)*100.0
    print("precip")
    precip=get_timeseries(files,"RAINNC",y,x)/3600.0
    precip[1:]=precip[1:]-precip[:-1]
    print("sw")
    sw=get_timeseries(files,"SWDOWN",y,x)
    print("lw")
    lw=get_timeseries(files,"GLW",y,x)
    
    print("times")
    times=read_times(filesearch,"Times")

    return Bunch(wind=windspeed,winddir=windspeed*0, t=t,rh=rh,p=p,precip=precip,sw=sw,lw=lw,times=times)

def read_forcing_step(filename,lastforcing=None):
    """Load a time series of WRF forcing data for the lat/lon position given"""
    
    u=swim_io.read_nc(filename,"U10").data
    v=swim_io.read_nc(filename,"V10").data
    windspeed=np.sqrt(u**2+v**2)
    t=swim_io.read_nc(filename,"T2").data
    p=swim_io.read_nc(filename,"PSFC").data/100.0
    q=swim_io.read_nc(filename,"Q2").data
    rh=units.mr2rh(t,p,q)*100.0
    precip=swim_io.read_nc(filename,"RAINNC").data/3600.0
    rawprecip=precip.copy()
    if lastforcing:
        precip-=lastforcing.rawprecip[-1,:,:]
    precip[1:,:,:]=precip[1:,:,:]-precip[:-1,:,:]
    sw=swim_io.read_nc(filename,"SWDOWN").data
    lw=swim_io.read_nc(filename,"GLW").data
    times=read_times(filename,"Times")

    return Bunch(wind=windspeed,winddir=windspeed*0, t=t,rh=rh,p=p,precip=precip,
                sw=sw,lw=lw,times=times,rawprecip=rawprecip)


    
def write_forcing(filename,forcing):
    outfile=open(filename,"au")
    
#   UTC date/time        windspeed       wind dir         temperature      humidity        pressure           shortwave      longwave       precipitation
#  1998 01 01 06 30     5.6300001144   178.0000000000   263.9499816895    86.0999984741  1002.0000000000     0.0000000000   281.0000000000     0.0000000000
    formatstring="{0.times[_i_]}  {0.wind[_i_]:15.9} {0.winddir[_i_]:15.9} {0.t[_i_]:15.9} {0.rh[_i_]:15.9} "\
                    "{0.p[_i_]:15.9} {0.sw[_i_]:15.9} {0.lw[_i_]:15.9} {0.precip[_i_]:15.9}\n"
    for i in range(forcing.wind.size):
        curformat=formatstring.replace("_i_",str(i))
        outfile.write(curformat.format(forcing))
    outfile.close

def append_forcing(filename,forcing,y,x):
    outfile=open(filename,"au")
    
#   UTC date/time        windspeed       wind dir         temperature      humidity        pressure           shortwave      longwave       precipitation
#  1998 01 01 06 30     5.6300001144   178.0000000000   263.9499816895    86.0999984741  1002.0000000000     0.0000000000   281.0000000000     0.0000000000
    xystring="["+str(y)+"]["+str(x)+"]"
    formatstring="{0.times[_i_]}  {0.wind[_i_]_xy_:15.9} {0.winddir[_i_]_xy_:15.9} {0.t[_i_]_xy_:15.9} {0.rh[_i_]_xy_:15.9} "\
                    "{0.p[_i_]_xy_:15.9} {0.sw[_i_]_xy_:15.9} {0.lw[_i_]_xy_:15.9} {0.precip[_i_]_xy_:15.9}\n"
    formatstring=formatstring.replace("_xy_",xystring)
    
    for i in range(forcing.wind.shape[0]):
        curformat=formatstring.replace("_i_",str(i))
        outfile.write(curformat.format(forcing))
    outfile.close
    

def write_outputfile(template,outputfile,parameters,forcing):
    """docstring for write_outputfile"""
    shutil.copyfile(template,outputfile)
    for s in searchstrings:
        replace_in_file(outputfile,s,str(parameters[s[2:-2]]))
    write_forcing(outputfile,forcing)

def setup_base_file(template,outputfile,parameters):
    """docstring for write_outputfile"""
    shutil.copyfile(template,outputfile)
    for s in searchstrings:
        replace_in_file(outputfile,s,str(parameters[s[2:-2]]))
    

def run_for_lat_lon(lat,lon):
    """Write a noahmp Simple Driver input file from WRF geogrid and outputfiles"""
    output_file="met_data_"+str(lat)+"_"+str(lon)+".dat"
    parameters=load_parameters(wrf_geo_grid_file,wrf_file_search,lat,lon)
    forcing=load_forcing_data(wrf_file_search,lat,lon)
    write_outputfile(template_file, output_file, parameters,forcing)

def run_for_grid():
    files=glob.glob(wrf_file_search)
    files.sort()
    files=files[318:418]
    basefile=files[0]
    xlat=swim_io.read_nc(basefile,"XLAT").data[0,...]
    xlon=swim_io.read_nc(basefile,"XLONG").data[0,...]
    xlat=xlat[30:210,100:220]
    xlon=xlon[30:210,100:220]
    outputfiles=[]
    positions=[]
    parameters=[]
    for lat,lon in zip(xlat.flat,xlon.flat):
        y,x=get_latlon_pos(basefile,lat,lon)
        print(y,x)
        positions.append((y,x))
        outputfiles.append("met_grid_data_"+str(lat)+"_"+str(lon)+".dat")
        # parameters.append(load_parameters(wrf_geo_grid_file,wrf_file_search,lat,lon,y,x))
        # setup_base_file(template_file,outputfiles[-1],parameters[-1])
        
    lastforcing=None
    for f in files:
        print(f)
        forcing=read_forcing_step(f,lastforcing)
        for pos,outf in zip(positions,outputfiles):
            append_forcing(outf,forcing,pos[0],pos[1])
        lastforcing=forcing

def main():
    run_for_grid()
    # run_for_lat_lon(40.0+2.0/60.0, -105.0-33.0/60.0)
    # 
    # basefile=glob.glob(wrf_file_search)[0]
    # xlat=swim_io.read_nc(basefile,"XLAT").data[0,...]
    # xlon=swim_io.read_nc(basefile,"XLONG").data[0,...]
    # for lat,lon in zip(xlat.flat,xlon.flat):
    #     run_for_lat_lon(lat,lon)
    #     parameters.append(load_parameters(wrf_geo_grid_file,wrf_file_search,lat,lon))
    # 

## TESTING:
# import wrf2noahmp
# reload(wrf2noahmp)
# from wrf2noahmp import *
# dates=calc_times(wrf_file_search)
# write_outputfile(template_file, output_file, parameters,forcing)
# forcing=load_forcing_data(wrf_file_search,lat,lon)
# parameters=load_parameters(wrf_geo_grid_file,wrf_file_search,lat,lon)
if __name__ == '__main__':
    main()