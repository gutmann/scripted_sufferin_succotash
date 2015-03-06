import numpy as np
import units

def interp_nearest(q,zin,zout):
    """docstring for interp_nearest"""
    qout=np.zeros(zout.shape)
    nz,ny,nx=zout.shape
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                i=0
                while (zin[i,y,x]<zout[z,y,x]):
                    i+=1
                    if i==zin.shape[0]:
                        i-=1
                        break
                i=max(i-1,0)
                qout[z,y,x]=q[i,y,x]
    return qout

def interp_pressure(pin,zin,zout,tin):
    """docstring for interp_pressure"""
    
    if np.any(np.array(zin.shape) != np.array(zout.shape)):
        z = interp_nearest(zin,zin,zout)
        p = interp_nearest(pin,zin,zout)
        t = interp_nearest(tin,zin,zout)
    else:
        z=zin
        p=pin
        t=tin
    
    outputp=np.zeros(p.shape)
    dtdz=np.zeros(t.shape)
    dtdz[:-1]=(t[1:]-t[:-1])/(z[1:]-z[:-1])
    dtdz[-1] = dtdz[-2]
    outputp=units.zt2p(zout-z,p0=p,t0=t,dtdz= -0.0065)
    return outputp

def _vinterp(data,zin,zout):
    nz,ny,nx=zout.shape
    nzin=zin.shape[0]
    output=np.zeros(zout.shape)
    # this should really be done as inline C
    for y in range(ny):
        for x in range(nx):
            i=0
            for z in range(nz):
                
                while (zout[z,y,x]<zin[i,y,x]):
                    i-=1
                    if i<0:
                        i=0
                        break
                
                while (zout[z,y,x]>zin[i,y,x]):
                    i=i+1
                    if i==nzin:
                        break
                
                if i==0:
                    output[z,y,x] = data[0,y,x]
                elif i==nzin:
                    output[z,y,x] = data[-1,y,x]
                else:
                    weight= (zout[z,y,x] - zin[i-1,y,x]) / (zin[i,y,x] - zin[i-1,y,x])
                    output[z,y,x] = weight * data[i,y,x] + (1-weight) * data[i-1,y,x]
    return output
                    

def interp(inputvar,inputz,outputz,vartype="data",inputt=None,inputp=None):
    """Interpolate an input 3D dataset from an input 3D elevation grid to the output 3D grid
    
    if vartype is set to "p" pressure, it will perform a more sophisticated interpolation, 
        but that requires REAL temperature data
    if vartype is set to "t" temperature, it will convert to potential temperature before
        converting, but that requires pressure data
    """
    # for pressure interpolation, use a temperature dependant vertical shift
    if vartype=="p":
        if inputt==None:
            raise(ValueError("Need Temperature to interpolate pressure"))
        return interp_pressure(inputvar,inputz,outputz,inputt)

    # if we were given real temperature, convert to potential temperature
    if vartype=="t":
        if inputp==None:
            raise(ValueError("Need Pressure data to interpolate temperature"))
            
        else:
            inputvar=inputvar*units.exner(inputp)
    

    output=_vinterp(inputvar,inputz,outputz)
            
            
            
    if vartype=="t":
        inputvar=inputvar/units.exner(inputp)
        print("WARNING!!  Returning potential temperature, not REAL temperature")
        print("  Use interpolated pressure to convert back to REAL temperature")
    
    return output


def _vinterp_multivar(data,output,timestep,zin,zout,varlist):
    # perform vertical interpolation as in _vinterp above, but across all variables
    # this way we don't have to find and compute weights multiple times
    
    nt,nz,ny,nx=data[varlist[0]].shape
    # this should really be done as inline C...
    # but I'm not sure that is easy with input/output classes like dictionaries (or even with a real dictionary)
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                # initial guess starts at the current level
                i=z
                # first make sure our initial guess is below the required level
                while (zout[z,y,x]<zin[i,y,x]):
                    i-=1
                    # make sure i>=0
                    if i<0:
                        i=0
                        break
                # now search upwards until we find the first level that is above the required level
                while (zout[z,y,x]>zin[i,y,x]):
                    i=i+1
                    # make sure i<=nz
                    if i==nz:
                        break
                # finally compute the weight based on the distance to the grid level below
                if (i!=0) and (i!=nz):
                    weight= (zout[z,y,x] - zin[i-1,y,x]) / (zin[i,y,x] - zin[i-1,y,x])
                # now loop through variables interpolating as appropriate
                for v in varlist:
                    if i==0:
                        output[v][timestep,z,y,x] = data[v][timestep,0,y,x]
                    elif i==nz:
                        output[v][timestep,z,y,x] = data[v][timestep,-1,y,x]
                    else:
                        output[v][timestep,z,y,x] = weight * data[v][timestep,i,y,x] + (1-weight) * data[v][timestep,i-1,y,x]


def interp_multivar(inputvar,outputvar,timestep,inputz,outputz,vartype="data",inputt=None,inputp=None):
    """Interpolate an input 3D dataset from an input 3D elevation grid to the output 3D grid
    
    if vartype is set to "p" pressure, it will perform a more sophisticated interpolation, 
        but that requires REAL temperature data
    if vartype is set to "t" temperature, it will convert to potential temperature before
        converting, but that requires pressure data
    """
    interp_keys=[]
    for vartype in inputvar.keys():
        # for pressure interpolation, use a temperature dependant vertical shift
        if vartype=="p":
            if inputt==None:
                raise(ValueError("Need Temperature to interpolate pressure"))
            outputvar[vartype][timestep]=interp_pressure(inputvar[vartype][timestep],
                                                         inputz, outputz, inputt)
        elif vartype!="z":
            interp_keys.append(vartype)

        # if we were given real temperature, convert to potential temperature
        if vartype=="t":
            if inputp==None:
                raise(ValueError("Need Pressure data to interpolate temperature"))
            
            else:
                inputvar[vartype][timestep]=inputvar[vartype][timestep]*units.exner(inputp)
                
    _vinterp_multivar(inputvar,outputvar,timestep,inputz,outputz,interp_keys)
    
    if "t" in inputvar.keys():
        inputvar["t"][timestep]=inputvar["t"][timestep]/units.exner(inputp)
        print("WARNING!!  Returning potential temperature, not REAL temperature")
        print("  Use interpolated pressure to convert back to REAL temperature")
    
