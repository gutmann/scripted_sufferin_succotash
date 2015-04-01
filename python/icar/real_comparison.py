#!/usr/bin/env python
import sys
import datetime
import glob
import operator as op

import numpy as np
import matplotlib.pyplot as plt
from bunch import Bunch

# import swm.post as sw
import mygis

wrf_dir="/glade/u/home/gutmann/scratch/wrfoutput/4km/2007/"

DIM_2D_SHAPE=3
DIM_3D_SHAPE=4

def echo(fn):
    def wrapped(*v, **k):
        print(fn.__name__)
        return fn(*v, **k)
    return wrapped

def exner(th,p):
    Rd=287.058
    cp=1004.0
    p0=100000
    pii=(p/p0)**(Rd/cp)
    return th * pii

class DataReader(object):
    # only directly accessible public attributes
    files=None
    times_per_file=1
    # curpos, last_rain, and last_rain_pos are accessible via @properties
    _curpos=0
    _pos_in_file=0
    _curfile=0
    
    _last_rain=None
    _last_rain_pos=-1
    _lr_pos_in_file=-1
    _lr_curfile=0
    _rainvar="RAINNC"
    _testvar=None
    
    # _var_names=["QVAPOR",[op.add,"QCLOUD","QICE"],"RAINNC",[op.add,"T2",300],"U","V","W"]
    _var_names=["QVAPOR",[op.add,"QCLOUD","QICE"],[op.add,"QRAIN","QSNOW"],[op.add,"T2",300],"U","V"]
    _short_names=dict(QVAPOR="qv",QCLOUD="qc",QICE="qc",RAINNC="rain",T="t",T2="t",U="u",V="v",W="w",QRAIN="rain",
                      qv="qv",qc="qc",rain="rain",qr="rain",th="t",u="u",v="v",w="w")
    _collapse_functions=dict(QVAPOR=np.mean,QCLOUD=np.sum,T=np.mean,U=np.mean,V=np.mean,W=np.mean,
                             QICE=np.sum,QRAIN=np.sum,QSNOW=np.sum,
                             qv=np.mean,qc=np.sum,th=np.mean,u=np.mean,v=np.mean,w=np.mean,
                             qi=np.sum,qr=np.sum,qs=np.sum,p=np.mean)

    _wrf_var_names=["QVAPOR",[op.add,"QCLOUD","QICE"],[op.add,"QRAIN","QSNOW"],"T2","U","V"]#[op.add,"T",290],"U","V"]
    _swm_var_names=["qv",[op.add,"qc","qi"],[op.add,"qr","qs"],[exner,"th","p"],"u","v"]

    
    x=slice(0,None) #by default take all data in the file in x,y, and z
    y=slice(0,None)
    # z=slice(0,None)
    z=slice(0,10)
    # zslices=dict(qv=slice(0,10),qc=slice(0,10),t=slice(1),)
    # yslices=dict()
    # yslices.setdefault(y)
    
    def __init__(self, filenames,start_pos=0,datatype="WRF"):
        super(DataReader,self).__init__()
        self.files=filenames
        
        self._datamodel=datatype
        if datatype=="WRF":
            self._var_names=self._wrf_var_names
            test_var=mygis.read_nc(self.files[0],self._var_names[0],returnNCvar=True)
            self.times_per_file=test_var.data.shape[0]
            test_var.ncfile.close()
            self.zaxis=0
            self.DIM_2D_SHAPE=3
            self.DIM_3D_SHAPE=4
        
        if datatype=="SWM":
            self._var_names=self._swm_var_names
            self.times_per_file=1
            self._rainvar="rain"
            tmp=self.y
            self.y=self.z
            self.z=tmp
            self.zaxis=1
            self.DIM_2D_SHAPE=2
            self.DIM_3D_SHAPE=3
        
        
        #note this calls the setter which will set pos_in_file and cur_file
        self.curpos=start_pos
    
    def _get_collapsing_func(self,varname):
        """docstring for get_collapsing_func"""
        try:
            myfunc=self._collapse_functions[varname]
        except:
            myfunc=np.mean
            
        return myfunc
    
    def collapse_z(self,data,varname):
        if len(data.shape)==3:
            myfunc=self._get_collapsing_func(varname)
            return myfunc(data,axis=self.zaxis)
        else:
            return data


    # Get/Set the position in the timeseries, while properly updating the filenumber and position in file
    @property
    def curpos(self):
        return self._curpos
        
    @curpos.setter
    def curpos(self,pos):
        self._curpos=pos
        self._pos_in_file= int(self._curpos) % int(self.times_per_file)
        self._curfile    = int(self._curpos) / int(self.times_per_file)


    # Get/Set the position in the timeseries, while properly updating the filenumber and position in file
    @property
    def last_rain_pos(self):
        return self._last_rain_pos
        
    @curpos.setter
    def last_rain_pos(self,pos):
        self._last_rain_pos=pos
        self._lr_pos_in_file= int(self._last_rain_pos) % int(self.times_per_file)
        self._lr_curfile    = int(self._last_rain_pos) / int(self.times_per_file)

                
    # Get/Set the last_rain variable
    @property
    def last_rain(self):
        if self._last_rain==None:
            self.last_rain_pos=self.curpos-1
            
            if (self._pos_in_file>0):
                nc_data=mygis.read_nc(self.files[self._curfile],self._rainvar,returnNCvar=True)
                self._last_rain=nc_data.data[self._last_rain_pos,self.y,self.x]
                nc_data.ncfile.close()
            elif (self._curfile==0):
                nc_data=mygis.read_nc(self.files[self._curfile],self._rainvar,returnNCvar=True)
                nx=nc_data.data.shape[1]
                ny=nc_data.data.shape[2]
                self._last_rain=np.zeros((nx,ny))[self.x,self.y]
                nc_data.ncfile.close()
            else:
                nc_data=mygis.read_nc(self.files[self._curfile-1],self._rainvar,returnNCvar=True)
                self._last_rain=nc_data.data[-1,self.x,self.y]
                nc_data.ncfile.close()
                
        # else: we already have a valid _last_rain, just return it this should be the case most of the time
        return self._last_rain
                
    
    @last_rain.setter
    def last_rain(self,value):
        if hasattr(value,__iter__):
            self.last_rain_pos=value[0]
            self._last_rain=value[1] 
        else:
            self.last_rain_pos=value
            self._last_rain=None # the getter will automagically generate last_rain
            
    
    def load_data(self,varname, filename=None, curtime=None):
        if type(varname)!=str:
            return varname
        
        if filename==None:
            filename=self.files[self._curfile]
        if curtime==None:
            curtime=self._pos_in_file
        
        data=mygis.read_nc(filename,varname,returnNCvar=True)
        dimlen=len(data.data.shape)
        # 2D vars e.g. RAINNC, rain
        if dimlen==self.DIM_2D_SHAPE:
            if dimlen==2:
                outputdata=data.data[self.y,self.x]
            else:
                outputdata=data.data[curtime,self.y,self.x]
        # 3D vars e.g. QVAPOR, qv
        elif dimlen==self.DIM_3D_SHAPE:
            if dimlen==3:
                outputdata=self.collapse_z(data.data[self.z,self.y,self.x],varname)
            else:
                outputdata=self.collapse_z(data.data[curtime,self.z,self.y,self.x],varname)
        else:
            raise IndexError("Do not know how to process {} dimensions".format(len(data.data.shape)))
            
        if varname==self._rainvar:
            curent_rain=outputdata[:]
            outputdata-=self.last_rain
            self.last_rain=(self.curpos,curent_rain)
        
        return outputdata
        
    def get_current_date(self):
        """Assumes a hard coded filename (e.g. WRF output filenames wrfout_d01_2007-01-01_00:00:00)"""
        if self._datamodel=="WRF":
            datestring=self.files[self._curfile].split("_")[2]+"-"+str(self._pos_in_file)
            return datetime.datetime.strptime(datestring,"%Y-%m-%d-%H")
        else:
            return datetime.datetime(2007,01,01,00)+datetime.timedelta(self.curpos/24.0)
    
    def __len__(self):
        return len(self.files)*self.times_per_file

    def __iter__(self):
        return self
        
    def __next__(self):
        
        self.curpos+=1
        output_data=Bunch()
        
        filename=self.files[self._curfile]
        for v in self._var_names:
            if type(v)==str:
                curdata=self.load_data(v)
                curvarname=v
            
            elif type(v)==list:
                cur_operator=v[0]
                for varname in v[1:]: 
                    if type(varname)==str:
                        curvarname=v[1]
                        break
                
                curdata=self.load_data(v[1])
                for curv in v[2:]:
                    next_data=self.load_data(curv)
                    cur_operator(curdata,next_data)
            
            output_data[self._short_names[curvarname]]=curdata
        
        output_data.date=self.get_current_date()
        return output_data
        
    next=__next__

clims=dict( qv=(0,0.004),
            qc=(0,0.0003),
            t=(260,310),
            u=(-15,15),
            v=(-15,15),
            rain=(0,0.000005))
def make_subplot(data,ny,nx,curplot,v,extra_title):
    plt.subplot(ny,nx,curplot)
    plt.imshow(data)
    plt.clim(clims[v])
    plt.colorbar()
    plt.title(v+extra_title)
    

def make_plots(data1,data2,date,fig=None):
    plt.close("all")
    if fig==None:
        fig=plt.figure(figsize=(24,14));
    else:
        fig.clear()
    ny=3
    nx=4
    curplot=0
    varnames=["qv","qc","u","v","t","rain"]
    for v in varnames:
        curplot+=1
        make_subplot(data1[v],ny,nx,curplot,v," "+str(date)[:14])
        curplot+=1
        make_subplot(data2[v],ny,nx,curplot,v," "+str(date)[:14])
    
    return fig
    

def main(swm_dir="output/",output_dir="./"):
    output_filename=output_dir+"vis_{}.png"
    wrf_files=glob.glob(wrf_dir+"wrfout*")
    wrf_files.sort()
    swm_files=glob.glob(swm_dir+"swim_out*")
    swm_files.sort()
    
    wrf_data=DataReader(wrf_files,datatype="WRF")
    swm_data=DataReader(swm_files,datatype="SWM")
    fig=plt.figure(figsize=(24,14));
    for i in range(len(wrf_data)):
        wrf=wrf_data.next()
        swm=swm_data.next()
        print(str(wrf.date),str(swm.date))
        sys.stdout.flush()
        
        fig=make_plots(swm,wrf,wrf.date,fig=fig)
        fig.savefig(output_filename.format(str(wrf.date).replace(" ","_")))
    


if __name__ == '__main__':

    global wrf_dir
    out_dir="./"
    swm_dir="output/"

    if len(sys.argv)>1:
        if sys.argv[1][:2]=="-h":
            print("Usage: real_comparison.py [swm_output_directory] [vis_output_directory] [wrf_dir]")
            sys.exit()
        
        swm_dir=sys.argv[1]
        if len(sys.argv)>2:
            out_dir=sys.argv[2]
        if len(sys.argv)>3:
            wrf_dir=sys.argv[3]

    
    main(swm_dir,out_dir)