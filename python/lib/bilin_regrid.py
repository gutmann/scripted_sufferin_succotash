'''The main external function is regrid, at the bottom'''
import numpy as np

def bilin_weights(yi,y,xi,x):
    '''Compute the bilinear weights for points surrounding xi,yi'''
    x0=np.abs((xi-x[0])/(x[1]-x[0]))
    x1=1-x0 # equivalent to np.abs((xi-x[1])/(x[1]-x[0]))
    x2=np.abs((xi-x[2])/(x[3]-x[2]))
    x3=1-x2 # equivalent to np.abs((xi-x[3])/(x[3]-x[2]))
    y5=y[0]*x1+y[1]*x0
    y6=y[2]*x3+y[3]*x2
    f1=(yi-y5)/(y6-y5)
    f2=1-f1# equivalent to (y6-yi)/(y6-y5)
    return np.array([x1*f2,x0*f2,x3*f1,x2*f1])

def load_geoLUT(lat1,lon1,lat2,lon2,subset=None):
    '''Create a Geographic Look Up Table
    
    lat/lon inputs should be grids of latitude and longitude
        lat1/lon1 = low resolution input grid
        lat2/lon2 = high resolution output grid
    '''
    # (lon2,lat2)=read_geo_latlon(georef,subset=subset)
    N=lat2.shape
    Nhi=lat1.shape
    # output data
    geoLUT=np.empty((N[0],N[1],4,3))

    # intermediate variables
    x=np.zeros(4).astype('i')
    y=np.zeros(4).astype('i')
    winhalfsize=5 # search window for next data point
    
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
                else:
                    x[1]=max(min(newx-dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx-dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy+dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy+dyinc,winsz[0]-1),0)
            else:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=max(min(newx+dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx+dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy-dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy-dyinc,winsz[0]-1),0)
                else:
                    x[1]=max(min(newx-dxinc,winsz[1]-1),0)
                    x[2]=newx
                    x[3]=max(min(newx-dxinc,winsz[1]-1),0)
                    y[1]=newy
                    y[2]=max(min(newy-dyinc,winsz[0]-1),0)
                    y[3]=max(min(newy-dyinc,winsz[0]-1),0)
            # finally compute the weights for each of the four surrounding points for a bilinear interpolation
            weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x])
            # store the current results in the output array
            geoLUT[j,i,:,0]=y+y0
            geoLUT[j,i,:,1]=x+x0
            geoLUT[j,i,:,2]=weights
    return geoLUT
    
def regrid(data,lat1,lon1,lat2,lon2,geoLUT=None,ymin=0,xmin=0):
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

    # if we weren't given a Geographic Look up table, create it now
    if not geoLUT:
        # if lat lon are linear instead of gridded, make them into grids
        if len(lat1.shape)==1:
            (lon1,lat1)=np.meshgrid(lon1,lat1)
        if len(lat2.shape)==1:
            (lon2,lat2)=np.meshgrid(lon2,lat2)
        geoLUT=load_geoLUT(lat1,lon1,lat2,lon2)
    
    # set up the output dataset 
    # WARNING: this could be HUGE, you may want to process the data in chunks so it will fit in memory
    outputdata=np.zeros((data.shape[0],geoLUT.shape[0],geoLUT.shape[1]),dtype=np.float32)
    # each of these four iterations corresponds to one of the four surrounding points
    # take that point, multiply by its weight and sum for bilinear interpolation
    for i in range(4):
        y=geoLUT[:,:,i,0].astype('i')-ymin
        x=geoLUT[:,:,i,1].astype('i')-xmin
        outputdata+=np.float32(data[:,y,x]*geoLUT[np.newaxis,:,:,i,2])
    return outputdata
