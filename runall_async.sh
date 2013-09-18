#!/bin/bash
# async_narr.py pr conus 12km NCEP &>SARoutput_pr_conus_12km &
# async_narr.py pr conus 6km NCEP &>SARoutput_pr_conus_6km &
# async_narr.py tasmin conus 12km NCEP &>SARoutput_tasmin_conus_12km &
# async_narr.py tasmin conus 6km NCEP &>SARoutput_tasmin_conus_6km &
# async_narr.py tasmax conus 12km NCEP &>SARoutput_tasmax_conus_12km &
# async_narr.py tasmax conus 6km NCEP &>SARoutput_tasmax_conus_6km 

# async_narr.py pr     conus 12km NCEP &>SARoutput_pr_conus_12km_ncep &
# async_narr.py tasmin conus 12km NCEP &>SARoutput_tasmin_conus_12km_ncep &
# async_narr.py tasmax conus 12km NCEP &>SARoutput_tasmax_conus_12km_ncep 

async_narr.py pr     conus 12km NARR &>SARoutput_pr_conus_12km_narr &
async_narr.py tasmin conus 12km NARR &>SARoutput_tasmin_conus_12km_narr &
async_narr.py tasmax conus 12km NARR &>SARoutput_tasmax_conus_12km_narr 

async_narr.py pr     conus 6km NARR &>SARoutput_pr_conus_6km_narr &
# for ((i=1;i<13;i++)); do async_narr.py pr conus 6km NARR $i &>SARoutput_pr_conus_6km_narr ; done
async_narr.py tasmin conus 6km NARR &>SARoutput_tasmin_conus_6km_narr &
# for ((i=1;i<13;i++)); do async_narr.py tasmin conus 6km NARR $i &>SARoutput_tasmin_conus_6km_narr ; done
async_narr.py tasmax conus 6km NARR &>SARoutput_tasmax_conus_6km_narr 
# for ((i=1;i<13;i++)); do async_narr.py tasmax conus 6km NARR $i &>SARoutput_tasmax_conus_6km_narr ; done

# async_narr.py pr     conus 6km NCEP &>SARoutput_pr_conus_6km_ncep 
# for ((i=1;i<13;i++)); do async_narr.py tasmin conus 6km NCEP $i &>SARoutput_tasmin_conus_6km_ncep ; done

# async_narr.py tasmax conus 6km NCEP &>SARoutput_tasmax_conus_6km_ncep



# async_narr.py pr e0 12km NCEP
# # async_narr.py pr e0 12km NARR
# async_narr.py pr e0 6km NCEP
# # async_narr.py pr e0 6km NARR
# async_narr.py pr e1 12km NCEP
# async_narr.py pr e1 12km NARR
# # async_narr.py pr e1 6km NCEP
# # async_narr.py pr e1 6km NARR
# async_narr.py tasmax e0 12km NCEP
# async_narr.py tasmax e0 12km NARR
# # async_narr.py tasmax e0 6km NCEP
# # async_narr.py tasmax e0 6km NARR
# async_narr.py tasmax e1 12km NCEP
# async_narr.py tasmax e1 12km NARR
# # async_narr.py tasmax e1 6km NCEP
# # async_narr.py tasmax e1 6km NARR
# async_narr.py tasmin e0 12km NCEP
# async_narr.py tasmin e0 12km NARR
# # async_narr.py tasmin e0 6km NCEP
# # async_narr.py tasmin e0 6km NARR
# async_narr.py tasmin e1 12km NCEP
# async_narr.py tasmin e1 12km NARR
# # async_narr.py tasmin e1 6km NCEP
# # async_narr.py tasmin e1 6km NARR
