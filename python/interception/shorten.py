#!/usr/bin/env python
import mygis
print("Reading data")
d=mygis.read_nc("Full_Movie.nc",returnNCvar=True)
data=d.data[:400,:,:,:]
print("Writing data")
mygis.write("short_movie.nc",data)
print("Done")