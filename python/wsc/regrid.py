import numpy as np

def norm(lat,lon):
    if len(lat.shape)==1:
        lon,lat=np.meshgrid(lon,lat)
    return lat,lon

def agg_lut(lathi,lonhi,latlo,lonlo):
    outputlut=np.zeros((lathi.shape[0],lathi.shape[1],2))-9999
    nx=lathi.shape[1]
    ny=lathi.shape[0]
    maxdist=((latlo[0,0]-latlo[1,1])**2 + (lonlo[0,0]-lonlo[1,1])**2)/1.9
    
    
    for i in range(ny):
        dists=(latlo-lathi[i,0])**2+(lonlo-lonhi[i,0])**2
        y,x=np.unravel_index(dists.argmin(),dists.shape)
        # print(i,ny,y,x,lathi[i,0],lonhi[i,0],latlo[y,x],lonlo[y,x])
        for j in range(nx):
            xmax=min(nx,x+3)
            xmin=max(0,x-3)
            ymax=min(ny,y+3)
            ymin=max(0,y-3)
            
            windists=((latlo[ymin:ymax,xmin:xmax]-lathi[i,j])**2+
                      (lonlo[ymin:ymax,xmin:xmax]-lonhi[i,j])**2)
            yoff,xoff=np.unravel_index(windists.argmin(),windists.shape)
            x=xoff+xmin
            y=yoff+ymin
            if windists[yoff,xoff]<maxdist:
                outputlut[i,j,:]=[y,x]
                # print(y,x,windists[yoff,xoff],maxdist)
    return outputlut
            
def aggdata(lut,data):
    arg


    xmins=lut[:,:,1].max(axis=0)
    ymins=lut[:,:,0].max(axis=1)
    
    tmp=np.where(xmins>0)[0]
    firstx=tmp[0]
    lastx=tmp[-1]
    newnx=lastx-firstx+2
    
    tmp=np.where(ymins>0)[0]
    firsty=tmp[0]
    lasty=tmp[-1]
    newny=lasty-firsty+2

    # print(newny,newnx)
    outputdata=np.zeros((data.shape[0],newny,newnx))
    n=np.zeros((1,newny,newnx))
    for i in range(data.shape[1]):
        for j in range(data.shape[2]):
            if lut[i,j,0]>0:
                outputdata[:,lut[i,j,0]-firsty,lut[i,j,1]-firstx]+=data[:,lut[i,j,0],lut[i,j,1]]
                n[0,lut[i,j,0]-firsty,lut[i,j,1]-firstx]+=1
                
    n[n==0]=-1
    return outputdata/n
    

def agg(data1,lat,lon,geo_lut=None):
    lat,lon=norm(lat,lon)
    lat1,lon1=norm(data1.lat,data1.lon)
    
    if geo_lut==None:geo_lut=agg_lut(lat1,lon1,lat,lon)
    
    return (geo_lut, aggdata(geo_lut,data1.data))
