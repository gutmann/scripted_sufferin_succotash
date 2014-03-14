# from netCDF4 import Dataset
import time
import os
import numpy as np
import Nio
from bunch import Bunch
import glob

def read_geo(filename):
    latnames=["lat","latitude","XLAT"]
    lonnames=["lon","longitude","XLONG"]
    for lat,lon in zip(latnames,lonnames):
        try:
            latdat=read_nc(filename,lat).data
            londat=read_nc(filename,lon).data
        except Exception as e:
            pass
    
    if londat.max()>180:
        londat=londat-360
    if len(londat.shape)==1:
        londat,latdat=np.meshgrid(londat,latdat)
    
    return Bunch(lat=latdat,lon=londat)

def Dataset(filename,mode="r",format="nc"):
    return Nio.open_file(filename,mode=mode,format=format)

def read_files(pattern,var="data",returnNCvar=False,axis=None):
    if type(pattern)==list:
        files=pattern
    else:
        files=glob.glob(pattern)
    files.sort()
    d=[]
    for f in files:
        d.append(read_nc(f,var=var,returnNCvar=returnNCvar).data)
    if axis!=None:
        d=np.concatenate(d,axis=axis)
    return d
    

def read_nc(filename,var="data",proj=None,returnNCvar=False):
    '''read a netCDF file and return the specified variable

    output is a structure :
        data:raw data as an array
        proj:string representation of the projection information
        atts:data attribute dictionary (if any)
    if (returnNCvar==True) then the Nio file is note closed and the Nio 
        representation of the variable is returned instead of being read into 
        memory immediately.  
    '''
    d=Nio.open_file(filename, mode='r',format="nc")
    outputdata=None
    if var != None:
        data=d.variables[var]
        attributes=d.variables[var].__dict__
        if returnNCvar:
            outputdata=data
        else:
            if len(data.shape)==0:
                outputdata=data.get_value()
            else:
                outputdata=data[:]
    outputproj=None
    if proj!=None:
        projection=d.variables[proj]
        outputproj=str(projection)
    
    
    if returnNCvar:
        return Bunch(data=outputdata,proj=outputproj,ncfile=d,atts=attributes)
    d.close()
    return Bunch(data=outputdata,proj=outputproj,atts=attributes)


def _write1d(NCfile,data,varname="data",units=None,dtype='f',dims=('x',),attributes=None):
    nx=data.size
    NCfile.create_dimension(dims[0], nx)
    NCfile.create_variable(varname,dtype,dims)
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units
    if attributes:
        for k in attributes.keys():
            NCfile.variables[varname].__setattr__(k,attributes[k])

def _write2d(NCfile,data,varname="data",units=None,dtype='f',dims=('y','x'),attributes=None):
    (ny,nx)=data.shape
    if dims[0]=="time":ny=0
    NCfile.create_dimension(dims[1], nx)
    NCfile.create_dimension(dims[0], ny)
    NCfile.create_variable(varname,dtype,dims)
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units
    if attributes:
        for k in attributes.keys():
            NCfile.variables[varname].__setattr__(k,attributes[k])

def _write3d(NCfile,data,varname="data",units=None,dtype='f',dims=('z','y','x'),attributes=None):
    (nz,ny,nx)=data.shape
    if dims[0]=="time":nz=0
    NCfile.create_dimension(dims[2], nx)
    NCfile.create_dimension(dims[1], ny)
    NCfile.create_dimension(dims[0], nz)
    NCfile.create_variable(varname,dtype,dims)
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units
    if attributes:
        for k in attributes.keys():
            NCfile.variables[varname].__setattr__(k,attributes[k])

def _write4d(NCfile,data,varname="data",units=None,dtype='f',dims=('t','z','y','x'),attributes=None):
    (nt,nz,ny,nx)=data.shape
    if dims[0]=="time":nt=0
    NCfile.create_dimension(dims[3], nx)
    NCfile.create_dimension(dims[2], ny)
    NCfile.create_dimension(dims[1], nz)
    NCfile.create_dimension(dims[0], nt)
    NCfile.create_variable(varname,dtype,dims)
    NCfile.variables[varname][:]=data.astype(dtype)
    if units!=None:
        NCfile.variables[varname].units=units
    if attributes:
        for k in attributes.keys():
            NCfile.variables[varname].__setattr__(k,attributes[k])


def addvar(NCfile,data,varname,dims,dtype='f',attributes=None):
    for i,d in enumerate(dims):
        if not(d in NCfile.dimensions):
            NCfile.create_dimension(d,data.shape[i])
    
    NCfile.create_variable(varname,dtype,dims)
    NCfile.variables[varname][:]=data.astype(dtype)
    if attributes:
        for k in attributes.keys():
            NCfile.variables[varname].__setattr__(k,attributes[k])

def write(filename,data,dtype='f',varname="data",dims=None,units=None,attributes=None,
          lat=None,lon=None,extravars=None,history=""):
    """write a netcdf file 
    
    filename = name of output netcdf file (.nc will be appended automatically)
    data = data to write to file
    dtype = data type of data (see below)
    varname = name of variable to create in the netcdf file
    units = units for the data
    lat = a latitude variable to add (here for legacy reasons, use extravars)
    lon = ditto    
    attribues: bunch/dictionary with key/values pairs to be added as attributes
    extravars is a list of variables to add
        each variable is a bunch or other class with attribues:
            data = data to write
            name = name of variable
            dims = dimensions to use: i.e. [[[[t],z],y],x] matching dimensions in primary data
            dtype= data type of output variable
                'd': 64 bit float
                'f': 32 bit float
                'l': long
                'i': 32 bit integer
                'h': 16 bit integer
                'b': 8 bit integer
                'S1': character
            attribues: bunch/dictionary with key/values pairs to be added as attributes
    """
    history = 'Created : ' + time.ctime() +'\nusing simple io.write by:'+os.environ['USER']+"  "+history
    NCfile=Nio.open_file(filename,mode="w",format="nc",history=history)
    if len(data.shape)==1:
        if dims==None:
           dims=('x',)
        _write1d(NCfile,data,varname=varname,units=units,dtype=dtype,dims=dims,attributes=attributes)
    if len(data.shape)==2:
        if dims==None:
           dims=('y','x')
        _write2d(NCfile,data,varname=varname,units=units,dtype=dtype,dims=dims,attributes=attributes)
    if len(data.shape)==3:
        if dims==None:
           dims=('z','y','x')
        _write3d(NCfile,data,varname=varname,units=units,dtype=dtype,dims=dims,attributes=attributes)
    if len(data.shape)==4:
        if dims==None:
           dims=('t','z','y','x')
        _write4d(NCfile,data,varname=varname,units=units,dtype=dtype,dims=dims,attributes=attributes)
    
    if lat!=None:
        if len(lat.shape)>1:
            NCfile.create_variable("lat",'f',(dims[-2],dims[-1]))
        else:
            NCfile.create_variable("lat",'f',(dims[-2],))
        NCfile.variables["lat"][:]=lat.astype('f')
    if lon!=None:
        if len(lon.shape)>1:
            NCfile.create_variable("lon",'f',(dims[-2],dims[-1]))
        else:
            NCfile.create_variable("lon",'f',(dims[-1],))
        NCfile.variables["lon"][:]=lon.astype('f')
    
    if extravars:
        for e in extravars:
            addvar(NCfile,e.data,e.name,e.dims,e.dtype,e.attributes)
    
    NCfile.close()


class NC_writer(object):
    NCfile=None
    curVar=None
    nx=0
    ny=0
    nz=None
    def __init__(self, filename,nx,ny,nz=None,var=None,dtype='f'):
        history = 'Created : ' + time.ctime() + '\nby:'+os.environ['USER']+" using NC_writer Class"
        self.NCfile=Nio.open_file(filename,mode='w',format="nc",history=history)
        self.NCfile.create_dimension('time', 0)
        self.NCfile.create_dimension('lat', ny)
        self.ny=ny
        self.NCfile.create_dimension('lon', nx)
        self.nx=nx
        if nz:
            self.NCfile.create_dimension('level',nz)
            self.nz=nz
        self.NCfile.create_variable('time','l',('time',))
        if var: self.addVar(var,dtype=dtype)
            
    
    def addVar(self,varname,dtype='f'):
        if self.NCfile:
            if self.nz:
                self.NCfile.create_variable(varname,dtype,('time','level','lat','lon'))
            else:
                self.NCfile.create_variable(varname,dtype,('time','lat','lon'))
            self.curVar=varname
    
    def appendToVar(self,data,varname=None,date=None,pos=None,dtype='f'):
        if varname==None:varname=self.curVar
        var=self.NCfile.variables[varname]
        if pos:
            n=pos
        else:
            n=self.NCfile.dimensions['time']
            if n==None:
                n=0
        if self.nz:
            var[n,:,:,:]=data.astype(dtype)
        else:
            var[n,:,:]=data.astype(dtype)
        if date:self.NCfile.variables['time'][n]=long(date)
    
    def close(self):
        if self.NCfile:
            self.NCfile.close()
            self.NCfile=None
            
    def __del__(self):
        self.close()
        # if self.NCfile:
        #     self.NCfile.close()
        #     self.NCfile=None
        # super(NC_writer,self).__del__()
            
