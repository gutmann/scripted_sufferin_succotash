import numpy as np
from bunch import Bunch

def dataXn(data,nsamples):
    """aggregate an NxM data array into an nsamples x nsamples array

    output=dataXn(data,nsamples)

    NOTE: if nsamples<0 it is treated as a stride length instead"""
    ny,nx=data.shape
    if nsamples>0:
        ystep=np.floor(ny/float(nsamples))
        xstep=np.floor(nx/float(nsamples))
        xsamples=nsamples
        ysamples=nsamples
    elif nsamples<0:
        xstep=-nsamples
        ystep=-nsamples
        xsamples=np.floor(nx/float(xstep))
        ysamples=np.floor(ny/float(ystep))

    outputdata=np.zeros((ysamples,xsamples))
    ystart=0
    for i in range(ysamples):
        xstart=0
        for j in range(xsamples):
            curdata=data[ystart:ystart+ystep,xstart:xstart+xstep]
            outputdata[i,j]=np.mean(curdata[np.isfinite(curdata)])
            xstart+=xstep
        ystart+=ystep
    return outputdata

def find_point(x,y,xmap,ymap,guess=None):
    if guess==None:
        dists = (xmap-x)**2 + (ymap-y)**2
        xout,yout = np.unravel_index(np.argmin(dists),dists.shape)

    else:
        window = 3
        nx = xmap.shape[1]
        xs = max(0,guess[0]-window)
        xe = guess[0]+window
        if xe>=nx:xe=None

        ny = xmap.shape[0]
        ys = max(0,guess[1]-window)
        ye = guess[1]+window
        if ye>=ny:ye=None

        dists = (xmap[ys:ye,xs:xe]-x)**2 + (ymap[ys:ye,xs:xe]-y)**2
        yout,xout = np.unravel_index(np.argmin(dists),dists.shape)
        xout+=xs
        yout+=ys

    guess=(xout,yout)

    max_delta = np.sqrt((xmap[0,0]-xmap[1,1])**2 + (ymap[0,0]-ymap[1,1])**2)/1.5

    curdist = np.sqrt((xmap[yout,xout]-x)**2 + (ymap[yout,xout]-y)**2)
    if (curdist > max_delta):
        xout=-1
        yout=-1

    return xout, yout, guess


def geoLUT(geoin, geoout):
    nx = geoin.lat.shape[1]
    ny = geoin.lat.shape[0]
# <<<<<<< HEAD
#     outputlut=np.zeros((ny,nx,3))
#
#     guess = None
#     for i in range(ny):
#         for j in range(nx):
#             x,y,guess = find_point(geoin.lon[i,j], geoin.lat[i,j], geoout.lon, geoout.lat, guess=guess)
#             outputlut[i,j,0] = x
#             outputlut[i,j,1] = y
#             outputlut[i,j,2]+= 1
#
#     return outputlut
# =======
    # outputlut=np.zeros((ny,nx,3))
    nxout = geoout.lat.shape[1]
    nyout = geoout.lat.shape[0]
    outputlut=np.empty(shape=(nyout,nxout),dtype=list)

    for i in range(ny):
        # reset the guess to None at the beginning of each row so it doesn't start
        # searching from the end of the last row.
        guess=None
        for j in range(nx):
            x,y,guess = find_point(geoin.lon[i,j], geoin.lat[i,j], geoout.lon, geoout.lat, guess=guess)
            # print(x,y,guess)
            if x>=0:
                if outputlut[y,x]==None:
                    outputlut[y,x] = [[i,j]]
                else:
                    outputlut[y,x].append([i,j])
            # this variation is easier to look at the output for testing
            # outputlut[i,j,0] = x
            # outputlut[i,j,1] = y
            # if (x>=0):
            #     outputlut[i,j,2]+= 1
    # for any points that didn't find any data, fill with the nearest neighbor
    for i in range(nyout):
        for j in range(nxout):
            if outputlut[i,j]==None:
                x,y,guess = find_point(geoout.lon[i,j],geoout.lat[i,j],geoin.lon, geoin.lat)
                outputlut[i,j]=[[guess[0],guess[1]]]

    return outputlut


def mode(inputdata):
    """Return the mode of a sequence"""
    # from : http://stackoverflow.com/questions/10797819/finding-the-mode-of-a-list-in-python
    from collections import Counter
    data = Counter(inputdata)
    return data.most_common(1)[0][0]  # Returns the most frequently occurring item


def agg_data(data, lut, method=np.mean):
    nx = lut.shape[1]
    ny = lut.shape[0]
    outputdata = np.zeros((ny,nx))

    for j in range(ny):
        for i in range(nx):
            if len(lut[j,i])>0:
                collection=[]
                for pos in lut[j,i]:
                    collection.append(data[pos[0],pos[1]])
                outputdata[j,i] = method(collection)

    return outputdata


def data2huc(data,hucmask,huc,minvalue=-100,maxvalue=1E5):
    """Compute mean of data where hucmask==huc

        if data is a 3D array loop over the first index and return a time series
        if data is a 2D array return the mean value

        also eliminates values <min or >max values (default=-100,1E5)

        usage: output=data2huc(data,hucmask,huc,minvalue=-100,maxvalue=1E5)
    """
    if len(data.shape)>2:
        # if data is a 3D array, loop over the first index... not ideal
        outputdata=np.zeros(data.shape[0])
        for i in range(data.shape[0]):
            curdata=data[i,hucmask==huc]
            outputdata[i]=np.mean(curdata[(curdata>=0) & (curdata<1E5)])
        return outputdata
    # if data is just a 2D array, simply return a mean value
    else:
        curdata=data[hucmask==huc]
        return np.mean(curdata[(curdata>=0) & (curdata<1E5)])

def data2hucs(data,hucmask,minvalue=-100,maxvalue=1E5):
    """Compute the mean of data for each huc in hucmask

    output=data2hucs(data,hucmask,minvalue=0,maxvalue=1E5)
        data = np.array(dims= 2 or 3)
        hucmask=np.array(dims=2) matching last 2 dims of data
        min/max values = valid range of values in data

    returns    a list of either mean values (if data==2D)
            or a list of np.arrays of len data.shape[0] (if data==3D)
    """
    huclist=np.unique(hucmask)
    huclist=huclist[huclist>0]
    hucdata=[]

    for huc in huclist:
        hucdata.append(data2huc(data,hucmask,huc,minvalue=minvalue,maxvalue=maxvalue))

    return Bunch(data=hucdata,hucs=huclist)
