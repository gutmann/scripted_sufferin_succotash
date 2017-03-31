#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import sys

def calc_window(pos,n, reg_n, search_window):
    """docstring for calc_window"""
    start=pos-search_window
    if start<reg_n:start=reg_n
    
    stop=pos+search_window
    if stop+reg_n>=n: stop=n-reg_n
    
    return (start, stop)

def compute_match(frameA, frameB,method="covariance",fsum=None):
    """docstring for compute_match"""
    if method=="correlation":
        return np.corrcoef(frameA.flat[:], frameB.flat[:])[0,1]
    if method=="difference":
        return 0-(np.abs(frameA-frameB)).sum()
        return 0-(np.abs(frameA.astype('f')-frameB)).sum()
    if method=="covariance":
        fb=frameB.astype('f')
        sum2=fb.sum()
        sum3=(frameA*fb).sum()
        return sum3-(fsum*sum2)/frameA.size
    

def find_match(frame_1, frame_2, loc_1, reg_nx=20, reg_ny=20, search_window=10):
    """Find a patch in frame_2 that matches the patch in frame_1
    
    frame_1 and frame_2 should be identically sized image data
    loc_1 should be a two element indexable object (x, y)
    
    nx, ny can define a width and height for the region
    """
    
    best_value=(-1,-1)
    best_match=None
    
    nx=frame_1.shape[0]
    ny=frame_1.shape[1]
    
    output_location=loc_1[:]
    
    # calculate the window to search through
    x_start, x_stop = calc_window(loc_1[0], nx, reg_nx, search_window)
    y_start, y_stop = calc_window(loc_1[1], ny, reg_ny, search_window)
    
    frame_match=frame_1[loc_1[0]-reg_nx:loc_1[0]+reg_nx,
                        loc_1[1]-reg_ny:loc_1[1]+reg_ny,:].astype('f')
    # fsum=None
    fsum=np.sum(frame_match)
    
    for x in range(x_start, x_stop):
        for y in range(y_start, y_stop):
            
            test_match=frame_2[x-reg_nx:x+reg_nx,
                                y-reg_ny:y+reg_ny,:]
            
            if best_match==None:
                best_match=compute_match(frame_match, test_match,fsum=fsum)
            else:
                curmatch=compute_match(frame_match, test_match,fsum=fsum)
            
                if curmatch>best_match:
                    output_location[:]=x,y
                    best_match=curmatch
    
    return output_location

def vid2nc(filename="/Users/gutmann/Desktop/IMG_2303.m4v",resolution=(1920,1080,3),n=-1, outputfile="movie_data.nc"):
    """docstring for vid2nc"""
    import mygis
    import video_reader as vr
    vid=vr.Video_Reader(filename,resolution)
    
    outputdata=[]
    if n>0:
        for i in range(n):
            outputdata.append(vid.next()[np.newaxis,:,:,:])
    else:
        i=0
        for v in vid:
            print("\rFrame : {}".format(i),end="")
            sys.stdout.flush()
            outputdata.append(v[np.newaxis,:,:,:])
            i+=1
    
    outputdata=np.concatenate(outputdata,axis=0)
    mygis.write(outputfile,outputdata)

def main(filename="/Users/gutmann/Desktop/IMG_2303.m4v",
         resolution=(1080,1920,3),  # resolution=(2160, 3840, 3)
         startloc=[1100, 300],
         fixedloc=None,
         n=220,
         patch_size=20,
         show_video=0):
         
    """docstring for main"""
    print("Example video tracking")
    import video_reader as vr
    import matplotlib.pyplot as plt
    
    
    loclist=[startloc]
    if fixedloc!=None:
        fixedlist=[fixedloc]
    
    vid=vr.Video_Reader(filename, resolution)
    
    d=vid.next()
    for i in range(n):
        d2=vid.next()
        newloc=find_match(d, d2, loclist[-1], reg_nx=patch_size, reg_ny=patch_size)
        loclist.append(newloc)
        
        if fixedloc!=None:
            fixedloc=find_match(d, d2, fixedlist[-1], reg_nx=patch_size*2, reg_ny=patch_size*2)
            fixedlist.append(fixedloc)

        d=d2
        
        xs=newloc[1]
        ys=newloc[0]
        
        if show_video>0:
            if i%show_video==0:
                plt.clf()
                plt.imshow(d2[ys-200:ys+200,xs-200:xs+200,:].astype("uint8"))
                plt.title(i)
                plt.plot(200,200,'x',color="red",markersize=30, markeredgewidth=3)
                if fixedloc!=None:
                    plt.plot(200+(fixedloc[1]-xs),200+(fixedloc[0]-ys),'x',
                             color="black",markersize=30, markeredgewidth=3)

                plt.ylim(300,0)
                plt.xlim(0,350)
                plt.draw()
    
    loclist=np.array(loclist)
    plt.figure(figsize=(20,12))
    plt.subplot(1,2,1)
    plt.plot(loclist[:,0])
    plt.plot(loclist[:,0],'x')
    plt.subplot(1,2,2)
    plt.plot(loclist[:,1])
    plt.plot(loclist[:,1],'x')
    
    if fixedloc!=None:
        fixedlist=np.array(fixedlist)
        plt.figure()
        plt.plot(loclist[:,0]-fixedlist[:,0])
        plt.plot(loclist[:,0]-fixedlist[:,0],'x')
        plt.figure()
        plt.plot(loclist[:,1]-fixedlist[:,1])
        plt.plot(loclist[:,1]-fixedlist[:,1],'x')
        
    
    plt.show()

if __name__ == '__main__':
    main(filename="GOPR2574.MP4", resolution=(2160, 3840, 3), startloc=[1100,300], n=100, show_video=1)