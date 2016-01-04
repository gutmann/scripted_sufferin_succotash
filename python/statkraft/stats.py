from __future__ import absolute_import, print_function, division

import numpy as np

def mean(data):
    return data.mean(axis=0)

def wetfrac(data, threshold=0):
    wetdays=np.zeros(data.shape)
    wetdays[data>threshold]=1
    return wetdays.mean(axis=0)

def sort(data):
    return np.sort(data,axis=0)

def percentile(data, percentile=0.99, dsort=None):
    if dsort==None:dsort=sort(data)
    if percentile>1:percentile/=100.0

    ntimes=data.shape[0]
    return dsort[ntimes*percentile,:,:]
