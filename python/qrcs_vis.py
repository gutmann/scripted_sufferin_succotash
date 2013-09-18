#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import swim_io as io

def qrcs_vis():
    qr=io.read_files("swim_qr_2000-10*","qr")
    qs=io.read_files("swim_qs_2000-10*","qs")
    qc=io.read_files("swim_qc_2000-10*","qc")

    plt.figure();
    qrmax=0.0003
    qcmax=0.0003
    qsmax=0.0003
    for i in range(400,600):
        print(i)
        qrimg=qr[i][:,:,50,75:125]/qrmax
        qcimg=qc[i][:,:,50,75:125]/qcmax
        qsimg=qs[i][:,:,50,75:125]/qsmax
        
        qrimg[qrimg>1]=1
        qsimg[qsimg>1]=1
        qcimg[qcimg>1]=1
        
        
        outputimg=np.ones((3,20,50))
        outputimg[:2,...]-=qcimg
        outputimg[1:,...]-=qrimg
        outputimg[0,...] -=qsimg[0,...]
        outputimg[2,...] -=qsimg[0,...]
        
        outputimg[outputimg<0]=0
        outputimg=(outputimg*255).astype("uint8")
        
        plt.clf()
        plt.imshow(outputimg.transpose([1,2,0]))
        plt.draw()
        plt.savefig("qcsr_{:03}.png".format(i))
        
if __name__ == '__main__':
    qrcs_vis()