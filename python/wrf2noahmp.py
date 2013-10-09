#!/usr/bin/env python
import numpy as np

import swim_io as io

def replace_in_file(filename,searchstring,outputstring):
    """replace searchstring with outputstring in the file filename
    
    Read all lines in filename, replaceing searchstring with outputstring 
    and writing to a new temporary file as you go, then copy temp back to filename
    """
    with open(filename,"ru") as f:
        with open(filename+".temporary","wu") as output:
            for l in f:
                l=l.replace(searchstring,outputstring)
                output.write(l)
    os.rename(filename+".temporary",filename)

def load_parameters(geogrid_file,lat,lon):
    """Read noahmp parameters from a wrf geogrid (netcdf) file
    
    Reads: 
        GREENFRAC,SOILTYP,VEGTYP, 
    """
    pass

def main():
    """Write a noahmp Simple Driver input file from WRF geogrid and outputfiles"""
    parameters=load_parameters(wrf_geo_grid_file,lat,lon)
    forcing=load_forcing_data(wrf_file_search,lat,lon)
    write_outputfile(template_file, output_file, parameters,forcing)

if __name__ == '__main__':
    main()