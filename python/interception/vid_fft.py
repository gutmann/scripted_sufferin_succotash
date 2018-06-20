#!/usr/bin/env python
# from __future__ import print_function
import sys, os
import glob

import numpy as np
import matplotlib.pyplot as plt

import video_reader as vr
import mygis

global verbose
verbose=True

global parallel
parallel=False

global number_of_fft
number_of_fft = 30 #30

global fps
fps = 30.0


def read_vid(filename="/Users/gutmann/Desktop/IMG_2303.m4v",resolution=(1080,1920,3),n=-1):
    """docstring for read_vid"""
    vid=vr.Video_Reader(filename,resolution)

    if verbose:print("Reading:"+filename)

    outputdata=[]
    if n>0:
        for i in range(n):
            outputdata.append(vid.next()[np.newaxis,:,:,:])
    else:
        i=0
        for v in vid:
            if verbose:print("\r  Frame : {}".format(i),end="")
            sys.stdout.flush()
            outputdata.append(v[np.newaxis,:,:,:])
            i+=1

    if verbose:print("\nFinished reading:"+filename)

    outputdata=np.concatenate(outputdata,axis=0)
    return outputdata

def compute_fft(i,data_subset):

    for_fft_data = data_subset[:,:,i:]

    freq_space = np.fft.fft(for_fft_data,axis=2)

    frequencies = np.fft.fftfreq(for_fft_data.shape[2], d=1/fps) [ 10:60 ]

    maxf = frequencies[np.argmax(np.abs(freq_space[:,:,10:60].real) + np.abs(freq_space[:,:,-10:-60:-1].real),axis=2)]
    f_ampl = np.max(np.abs(freq_space[:,:,10:60].real) + np.abs(freq_space[:,:,-10:-60:-1].real),axis=2)

    return maxf, f_ampl

if parallel:
    from multiprocessing import Pool
    n_processors = 10
    p = Pool(n_processors)


def compute_frequency(data):
    if verbose:print("Computing frequencies")

    xmin=700
    xmax=1200
    # xmin=330
    # xmax=760
    # xmin=800
    # xmax=1050
    xstep = 1 if (xmax > xmin) else -1

    ymin=750
    ymax=450
    # ymin=1080
    # ymax=1080-400
    # ymax=1080-200
    ystep = 1 if (ymax > ymin) else -1

    data_subset = np.transpose(data[:,ymin:ymax:ystep,xmin:xmax:xstep,0],axes=[1,2,0])

    maxf = np.zeros((number_of_fft, data_subset.shape[0], data_subset.shape[1]))
    f_ampl = np.zeros((number_of_fft, data_subset.shape[0], data_subset.shape[1]))

    if parallel:
        res=[]

    for i in range(number_of_fft):

        if parallel:
            if verbose: print("\r Starting FFT:{:3}".format(i+1), end="")
            res.append(p.apply_async(compute_fft, (i,data_subset)))
        else:
            f,a = compute_fft(i,data_subset)
            if verbose: print("\r Completed FFT:{:3}".format(i+1), end="")
            maxf[i] = f
            f_ampl[i] = a
        sys.stdout.flush()

    if parallel:
        if verbose:print("  Getting FFT results")
        for i in range(number_of_fft):
            if verbose: print("\r Completed FFT:{:3}".format(i+1), end="")
            sys.stdout.flush()
            maxf[i], f_ampl[i] = res[i].get()

    if verbose:print("\r Completed all FFT computations")


    # tmp=np.argmax(f_ampl,axis=0).astype("uint8")
    # bestf=np.zeros(tmp.shape)
    # for i in range(tmp.shape[0]):
    #         for j in range(tmp.shape[1]):
    #                 bestf[i,j]=maxf[tmp[i,j],i,j]
    #
    # besta = f_ampl.max(axis=0)
    # bestf = np.ma.array(bestf, mask = besta<1500)
    # return (bestf, besta)

    return (maxf, f_ampl)

def write(filename, data):
    mygis.write(filename,data)

def make_image(data, inputfile):
    plt.clf()
    plt.imshow(data[0].astype('uint8'),origin="upper")
    plt.savefig(inputfile.split(".")[0]+".png")


    plt.clf()
    print(data.shape)
    plt.imshow(data[:,:,:,1].std(axis=0),origin="upper")
    plt.colorbar()
    plt.clim(5,25)
    plt.savefig(inputfile.split(".")[0]+"_std.png")


def complete(filename, outputdir=""):
    return os.path.isfile(outputdir+"ampl_"+filename.split("/")[-1][:-4]+".nc")


def main(filesearch, outputdir="output3/"):
    files=glob.glob(filesearch)
    files.sort()
    for f in files:
        if verbose:print(f)
        if not complete(f, outputdir):
            # print("Not Complete:"+f.split("/")[-1])
            data = read_vid(filename=f)
            make_image(data, f)
            freq, ampl = compute_frequency(data)

            write(outputdir+"freq_"+f.split("/")[-1][:-4]+".nc", freq)
            write(outputdir+"ampl_"+f.split("/")[-1][:-4]+".nc", ampl)
        else:
            print("Already completed:"+f.split("/")[-1])


if __name__ == '__main__':
    verbose=True
    # main("niwot_gopro_20160603/*.MP4")
    if len(sys.argv)>1:
        # print(sys.argv[1])
        main(sys.argv[1]+"/*.MP4", outputdir="output/")
    else:
        print("ERROR: Please specify input directory")
