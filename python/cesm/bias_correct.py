25

import numpy as np
import glob
import mygis
from bunch import Bunch
import units

# assumes CESM-LE climo file generated simply by ncea averaging all historical years and all ensemble members
# assumes era-interim climatology (from trude, maybe from kyoko originally) has been regridded horizontally append
#   vertically to match the CESM-LE climo file.  but CESM-LE climo file is 6hrly, ERA-i file is daily. 
#   To do this, put U and V in their own files using ncks and rename their XLAT_U, etc to XLAT (so read_geo will find)
#   then regrid all mass and wind variables horizontally using regrid_file2file(regrid_vertical=False)
#   Then recombine U and V variables with mass variables using ncks. 
#   Finally vertically interpolate new file to the cesm_le climo using regrid_file2file(regrid_horizontal=False)

cesm_le_climo = "cesm_le_climo.nc"
erai_climo    = "regrid_master_era.nc"

era_data = Bunch()
cesm_data= Bunch()

era_variables = Bunch(u="UU",v="VV",t="TT",p="PRES",rh="RH",tskin="SKINTEMP")
cesm_variables= Bunch(qv="qv",u="u",v="v",th="theta",p="p",sw="swdown",lw="lwdown",
                      dz="dz",z="z",terrain="hgt",lat="lat",lon="lon",landmask="xland")

filesearch = "icar_001_20TR_1990-01-01_00:00:00"


biases = None

def running_ave(data):
    window = 15 # plus or minus 15 days
    outputdata=np.zeros(data.shape)
    nz = data.shape[0]
    for i in range(data.shape[0]):
        if i<window:
            outputdata[i] = data[:i+window+1].sum(axis=0)
            outputdata[i]+= data[i-window:].sum(axis=0)
            outputdata[i]/= window*2+1
        elif (i+window+1)>nz:
            outputdata[i] = data[:(i+window+1)-nz].sum(axis=0)
            outputdata[i]+= data[i-window:].sum(axis=0)
            outputdata[i]/= window*2+1
        else:
            outputdata[i] = data[i-window:i+window+1].mean(axis=0)
    return outputdata

def cesm_running_ave(data):
    nx=data.shape[-1]
    ny=data.shape[-2]
    
    if len(data.shape)==3:
        reshaped_data = data.reshape((365,4,ny,nx))
        
    elif len(data.shape)==4:
        nz = data.shape[-3]
        reshaped_data = data.reshape((365,4,nz,ny,nx))
        
    return running_ave(reshaped_data.mean(axis=1) )

def write_biases(fileprefix):
    for k in biases.keys():
        mygis.write(fileprefix+k+".nc",biases[k])

def compute_bias():
    global biases
    
    biases=Bunch()

    easy_vars = ["u","v","p"] # eventually should add tskin
    for v in easy_vars:
        # compute bias in U, V, and pressure
        era_data[v]  = running_ave( mygis.read_nc(erai_climo,era_variables[v]).data )
        cesm_data[v] = cesm_running_ave( mygis.read_nc(cesm_le_climo,cesm_variables[v]).data )
        biases[v] = era_data[v] - cesm_data[v]

    # compute bias in potential temperature
    era_exner    = units.exner( mygis.read_nc(erai_climo,era_variables.p).data/100)
    era_data.th  = running_ave( mygis.read_nc(erai_climo,era_variables.t).data / era_exner )
    era_data.t   = running_ave( mygis.read_nc(erai_climo,era_variables.t).data )
    cesm_exner   = units.exner( mygis.read_nc(cesm_le_climo,cesm_variables.p).data/100)
    cesm_data.th = cesm_running_ave( mygis.read_nc(cesm_le_climo, cesm_variables.th).data )
    cesm_data.t  = cesm_running_ave( mygis.read_nc(cesm_le_climo, cesm_variables.th).data * cesm_exner )
    biases.th    = era_data.th - cesm_data.th
    
    # compute bias in relative humidity
    era_data.rh  = running_ave( mygis.read_nc(erai_climo,era_variables.rh).data)/100.0
    t = mygis.read_nc(cesm_le_climo, cesm_variables.th).data * cesm_exner
    p = mygis.read_nc(cesm_le_climo, cesm_variables.p).data / 100
    sh= mygis.read_nc(cesm_le_climo, cesm_variables.qv).data
    cesm_data.rh = cesm_running_ave( units.sh2rh(t,p,sh) )
    biases.rh = era_data.rh - cesm_data.rh

    print("Writing Biases")
    write_biases("bias_")

def write_bias_corrected_data(data, filename, global_atts=None):
    extra_vars=[]
    for k in data.keys():
        if k[-5:]!="_atts":
            if len(data[k].shape)==1:
                dims=("lev",)
            elif len(data[k].shape)==2:
                dims=("lat","lon")
            elif len(data[k].shape)==3:
                dims=("time","lat","lon")
            elif len(data[k].shape)==4:
                dims=("time","lev","lat","lon")
            print(dims)
            print(data[k].shape)
            extra_vars.append( Bunch(data=data[k],name=k,dims=dims,dtype="f",attributes=data[k+"_atts"] ) )
    
    d = extra_vars.pop()
    
    mygis.write(filename,d.data, varname=d.name,dims=d.dims, attributes=d.attributes, 
                 global_attributes=global_atts, extravars=extra_vars)

def main(filename):
    
    if biases==None:
        compute_bias()
    
    all_data = Bunch()
    print("Reading main data")
    for varname in cesm_variables.keys():
        all_data[varname] = mygis.read_nc(filename,cesm_variables[varname]).data
        all_data[varname+"_atts"] = mygis.read_atts(filename,cesm_variables[varname])
    
    rh = units.sh2rh( all_data.th * units.exner(all_data.p/100), all_data.p/100, all_data.qv )
    
    print("Correcting Biases")
    nbiases = biases.u.shape[0]
    for i in range(all_data.qv.shape[0]):
        # because input data are 6hrly and biases are daily we compute a linear interpolation
        # between bias correction dates
        first = int( np.floor( (i-2)/4.0 ))
        last  = int( np.floor( (i+2)/4.0 )) % nbiases
        last_weight= ((i+2.5)%4)/4.0
        first_weight = 1-last_weight
        
        all_data.u[i]  += biases.u[last]*last_weight  + biases.u[first]*first_weight
        all_data.v[i]  += biases.v[last]*last_weight  + biases.v[first]*first_weight
        all_data.p[i]  += biases.p[last]*last_weight  + biases.p[first]*first_weight
        all_data.th[i] += biases.th[last]*last_weight + biases.th[first]*first_weight
        rh[i] += biases.rh[last]*last_weight + biases.rh[first]*first_weight
        rh[i][rh[i] < 1e-10] = 1e-10
        exner = units.exner(all_data.p[i]/100)
        all_data.qv[i] = units.rh2sh(all_data.th[i] * exner, all_data.p[i], rh[i])
    
    global_atts = mygis.read_atts(filename, global_atts=True)
    print("Writing data")
    write_bias_corrected_data(all_data, "bc_"+filename, global_atts=global_atts)
        

        

if __name__ == '__main__':
    all_files = glob.glob(filesearch)
    all_files.sort()
    for f in all_files:
        print("Correcting: "+f)
        main(f)
