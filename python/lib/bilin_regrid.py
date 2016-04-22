'''The main external function is regrid, at the bottom'''
import numpy as np
from bunch import Bunch

def bilin_weights(yi,y,xi,x,mask=None,xj=None,yj=None,case=None):
    '''Compute the bilinear weights for points surrounding xi,yi'''
    x0=np.abs((xi-x[0])/(x[1]-x[0]))
    x1=1-x0 # equivalent to np.abs((xi-x[1])/(x[1]-x[0]))
    x2=np.abs((xi-x[2])/(x[3]-x[2]))
    x3=1-x2 # equivalent to np.abs((xi-x[3])/(x[3]-x[2]))
    y5=y[0]*x1+y[1]*x0
    y6=y[2]*x3+y[3]*x2
    if (y6-y5)==0:
        raise ValueError("Divide by zero, increase window size")
        # print(y)
        # print(x1,x0,x3,x2)
        # print(x)
        # print(yi,xi)
        # print(xj,yj)
        # print(case)
        # print("-"*20)
        
    f1=(yi-y5)/(y6-y5)
    f2=1-f1# equivalent to (y6-yi)/(y6-y5)
    if mask==None:
        return np.array([x1*f2,x0*f2,x3*f1,x2*f1])
    else:
        weights=np.array([x1*f2,x0*f2,x3*f1,x2*f1])
        weights[mask]=0
        if weights.max()>0:
            goodpoints=weights>0
            weights[goodpoints]/=np.sum(weights[goodpoints])
        return weights

def load_geoLUT(lat1,lon1,lat2,lon2,subset=None,mask=None,winhalfsize=5):
    '''Create a Geographic Look Up Table
    
    lat/lon inputs should be grids of latitude and longitude
        lat1/lon1 = low resolution input grid
        lat2/lon2 = high resolution output grid
    '''
    # winhalfsize=7 # search window for next data point 5 is faster, but 7 is safer - one case failed on 5
    # (lon2,lat2)=read_geo_latlon(georef,subset=subset)
    if len(lat1.shape)==1:
        lon1,lat1=np.meshgrid(lon1,lat1)
    if len(lat2.shape)==1:
        lon2,lat2=np.meshgrid(lon2,lat2)
    N=lat2.shape
    Nhi=lat1.shape
    # output data
    geoLUT=np.empty((N[0],N[1],4,3))

    # intermediate variables
    x=np.zeros(4).astype('i')
    y=np.zeros(4).astype('i')
    
    # figure out which direction is positive in latitude and longitude (just in case)
    dxinc=np.sign(lon1[1,1]-lon1[0,0]).astype('i')
    dyinc=np.sign(lat1[1,1]-lat1[0,0]).astype('i')

    # loop over y
    for i in range(N[1]):
        # find the first positions so the rest can be done faster relative to it
        j=0
        # compute distances for the ENTIRE low-res grid for the first input point
        dists=(lat1-lat2[j,i])**2 + (lon1-lon2[j,i])**2
        # find the minimum distance position
        (lasty,lastx)=np.unravel_index(dists.argmin(), dists.shape)
        (prevx,prevy)=(lastx,lasty)
        # create a window around that position to start searching as we loop through the rest of it
        y0=max(lasty-winhalfsize,0)
        y1=min(lasty+winhalfsize,Nhi[0])
        x0=max(lastx-winhalfsize,0)
        x1=min(lastx+winhalfsize,Nhi[1])
        latwin=lat1[y0:y1,x0:x1]
        lonwin=lon1[y0:y1,x0:x1]
        for j in range(N[0]):
            # if we have moved update the window position
            if (prevx!=lastx) or (prevy!=lasty):
                y0=max(lasty-winhalfsize,0)
                y1=min(lasty+winhalfsize,Nhi[0])
                x0=max(lastx-winhalfsize,0)
                x1=min(lastx+winhalfsize,Nhi[1])
                latwin=lat1[y0:y1,x0:x1]
                lonwin=lon1[y0:y1,x0:x1]
                (prevx,prevy)=(lastx,lasty)
            # compute distances just over this window
            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
            # find minimum distance position
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            # convert that position to a position in the large grid
            lastx=newx+x0
            lasty=newy+y0
            # this is the first point we care about
            x[0]=newx
            y[0]=newy
            # now find the position of the other three surrounding points
            winsz=latwin.shape
            if latwin[newy,newx]<lat2[j,i]:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=max(min(newx+dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx+dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy+dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy+dyinc,winsz[0]-1),0)
                    case=0
                else:
                    x[1]=max(min(newx-dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx-dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy+dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy+dyinc,winsz[0]-1),0)
                    case=1
            else:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=max(min(newx+dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx+dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy-dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy-dyinc,winsz[0]-1),0)
                    case=2
                else:
                    x[1]=max(min(newx-dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx-dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy-dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy-dyinc,winsz[0]-1),0)
                    case=3
            # finally compute the weights for each of the four surrounding points for a bilinear interpolation
            if mask!=None:
                weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x],mask=mask[y+y0,x+x0],xj=x,yj=y,case=case)
            else:
                weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x])
            # store the current results in the output array
            geoLUT[j,i,:,0]=y+y0
            geoLUT[j,i,:,1]=x+x0
            geoLUT[j,i,:,2]=weights
    return geoLUT
    
def regrid(data,lat1=None,lon1=None,lat2=None,lon2=None,geoLUT=None,ymin=0,xmin=0,output_geoLUT=False,missing=None):
    '''Regrid a dataset from the lat1,lon1 grid to the lat2,lon2 grid
    
    data: low resolution data to be regridded
            np.array(ntimes x nlat x nlon)
    lat/lon1: input low resolution grid 
            np.array(nlat x nlon) or np.array(nlat),np.array(nlon)
    lat/lon2: output high resolution grid
            np.array(nlat x nlon) or np.array(nlat),np.array(nlon)
    NOTE: lon1,lon2 must be in the same system (e.g. both 0:360 or both -180:180)
        and lat2,lon2 must be geographically internal to lat1,lon1 (may need some buffer too)
    WARNING: make sure that the output grid will fit in memory or this will take a LONG time
        if the output grid is too short, process it in N time slices.  
        Calculate geoLUT first and pass it on each iteration to save time
    '''
    if missing!=None:
        mask=data[0,...]==missing
    else:
        mask=None

    # if we weren't given a Geographic Look up table, create it now
    if geoLUT==None:
        # if lat lon are linear instead of gridded, make them into grids
        if len(lat1.shape)==1:
            (lon1,lat1)=np.meshgrid(lon1,lat1)
        if len(lat2.shape)==1:
            (lon2,lat2)=np.meshgrid(lon2,lat2)
        geoLUT=load_geoLUT(lat1,lon1,lat2,lon2,mask=mask)
        
    # set up the output dataset 
    # WARNING: this could be HUGE, you may want to process the data in chunks so it will fit in memory
    if len(data.shape)==2:
        data=data.reshape((1,data.shape[0],data.shape[1]))
    outputdata=np.zeros((data.shape[0],geoLUT.shape[0],geoLUT.shape[1]),dtype=np.float32)
    # each of these four iterations corresponds to one of the four surrounding points
    # take that point, multiply by its weight and sum for bilinear interpolation
    for i in range(4):
        y=geoLUT[:,:,i,0].astype('i')-ymin
        x=geoLUT[:,:,i,1].astype('i')-xmin
        outputdata+=np.float32(data[:,y,x]*geoLUT[np.newaxis,:,:,i,2])
    
    if output_geoLUT:
        return outputdata,geoLUT
    else:
        return outputdata

def get_geometry(z):
    """docstring for get_geometry"""
    if z[0]>z[-1]:
        return Bunch(bottom=len(z)-1,top=0,step=-1)
    else:
        return Bunch(bottom=0,top=len(z)-1,step=1)
        

def vLUT_1d_to_1d(inputz,outputz,geom_in=None,geom_out=None):
    """docstring for vLUT_1d_to_1d"""
    
    outputLUT=np.zeros((outputz.shape[0],3)) # 3 elements are point 1, point2, weight for point1 (w2=1-w1)
    
    if geom_in==None:
        geom_in=get_geometry(inputz)
    if geom_out==None:
        geom_out=get_geometry(outputz)
        
    for i in range(geom_out.bottom,geom_out.top+geom_out.step,geom_out.step):
        curz=outputz[i]
        if inputz[geom_in.bottom]>=curz:
            outputLUT[i,:]=np.array([geom_in.bottom,0,1.0])
        elif inputz[geom_in.top]<=curz:
            outputLUT[i,:]=np.array([geom_in.top,0,1.0])
        else:
            for z in range(geom_in.bottom,geom_in.top,geom_in.step):
                nextz=z+geom_in.step
                if inputz[nextz]>curz:
                    weight=1-(curz-inputz[z])/(inputz[nextz]-inputz[z])
                    outputLUT[i,:]=np.array([z,nextz,weight])
                    break
                    
    return outputLUT

def vLUT_1d_to_3d(inputz,outputz):
    outputshape=np.zeros(4)
    outputshape[:3]=outputz.shape
    outputshape[3]=3
    outputLUT=np.zeros(outputshape)
    
    for i in range(outputz.shape[1]):
        for j in range(outputz.shape[2]):
            outputLUT[:,i,j,:]=vLUT_1d_to_1d(inputz,outputz[:,i,j])
    
    return outputLUT

def vLUT_3d_to_3d(inputz,outputz):
    for i in range(1,3):
        if inputz.shape[i]!=outputz.shape[i]:
            raise ValueError("Shapes must match in last two dimensions.")
            
    outputshape=np.zeros(4)
    outputshape[:3]=outputz.shape
    outputshape[3]=3
    outputLUT=np.zeros(outputshape)
    
    for i in range(outputz.shape[1]):
        for j in range(outputz.shape[2]):
            outputLUT[:,i,j,:]=vLUT_1d_to_1d(inputz[:,i,j],outputz[:,i,j])
    
    return outputLUT
    
        
def load_vLUT(inputz,outputz):
    """docstring for load_vLUT"""
    if len(inputz.shape)==1:
        if len(outputz.shape)==3:
            return vLUT_1d_to_3d(inputz,outputz)
        if len(outputz.shape)==1:
            return vLUT_1d_to_1d(inputz,outputz)

    if len(inputz.shape)==3:
        if len(outputz.shape)==3:
            return vLUT_3d_to_3d(inputz,outputz)
    
    else:
        raise ValueError("Invalid dimensions (valid dimensionality:1D->1D, 1D->3D, 3D->3D)")


def vinterp(data,inputz=None,outputz=None,vLUT=None):
    """docstring for vinterp(era_on_gcm_grid,vLUT=pLUT)"""
    if vLUT==None:
        vLUT=load_vLUT(inputz,outputz)
        
    if len(vLUT.shape)==4:
        outputdata=np.zeros(vLUT.shape[:-1])
        nz,ny,nx=outputdata.shape
        for k in range(nz):
            for i in range(ny):
                for j in range(nx):
                    z1=vLUT[k,i,j,0]
                    z2=vLUT[k,i,j,1]
                    w=vLUT[k,i,j,2]
                    outputdata[k,i,j]=data[z1,i,j]*w+data[z2,i,j]*(1-w)
                    
    elif len(vLUT.shape)==2:
        if len(data.shape)==3:
            outputdata=np.zeros((vLUT.shape[0],data.shape[1],data.shape[2]))
            nz,ny,nx=outputdata.shape
            for k in range(nz):
                for i in range(ny):
                    for j in range(nx):
                        z1=vLUT[k,0]
                        z2=vLUT[k,1]
                        w=vLUT[k,2]
                        outputdata[k,i,j]=data[z1,i,j]*w+data[z2,i,j]*(1-w)
        if len(data.shape)==1:
            outputdata=np.zeros(vLUT.shape[0])
            nz=outputdata.shape
            for k in range(nz):
                z1=vLUT[k,0]
                z2=vLUT[k,1]
                w=vLUT[k,2]
                outputdata[k]=data[z1]*w+data[z2]*(1-w)
            
    else:
        raise ValueError("Invalid vLUT dimensions 1D or 3D only")
    
    
    return outputdata
    
    
