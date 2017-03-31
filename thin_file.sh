#!/bin/sh

# ncks -v rain,rain_rate,snow,crain,ta2m,hus2m,u10m,v10m,soil_w,soil_t,snw,rsds,rlds $1 short_$1
ncks -v time,lat,lon,rain,rain_rate,snow,graupel,crain,rsds,rlds,ts,ta2m,hus2m,u10m,v10m,hfss,hfls,rlus,hfgs,soil_w,soil_t,snw,canwat $1 short_$1
if [[ -e short_$1 ]]; then
    rm $1
fi
