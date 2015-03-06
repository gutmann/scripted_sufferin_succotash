#!/usr/bin/env python
# encoding: utf-8
"""
get_chimney_ndvi.py

Created by Ethan Gutmann on 2011-06-20.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""

import numpy as np
import mygis
import glob
from bunch import Bunch
import re

def read_gain_offset(filename):
	DNrange=255-1
	maxvals=np.zeros(9) #technically there are 9 potential gain/offset pair
	minvals=np.zeros(9)
	
	with open(filename,'r') as f:
		for line in f:
			if re.match(' +LM[AI?][XN?]_BAND[1-5] = .*',line):
				data=line.split(' = ')
				bandnum=np.int(data[0][-1])
				maxmin=data[0].split('_')[0][-1]
				if maxmin=='X':
					maxvals[bandnum-1]=np.float(data[1])
				else:
					minvals[bandnum-1]=np.float(data[1])
	gains=(maxvals-minvals)/DNrange
	offsets=minvals-gains
	return (gains,offsets)

def get_ndvi(utmpoint=None,niwot=False,chimney=True):
	b2file=glob.glob('*B20.TIF')[0]
	b3file=glob.glob('*B30.TIF')[0]
	b4file=glob.glob('*B40.TIF')[0]
	metafile=glob.glob('*MTL.txt')[0]
	Esun=np.array((1969,1840,1551,1041,225.7,82.07,1368))
	
	b2=mygis.read_tiff(b2file)
	b3=mygis.read_tiff(b3file)
	b4=mygis.read_tiff(b4file)
	
	(gains,offsets)=read_gain_offset(metafile)
	b2_dark_object=min(b2.data[np.where(b2.data>0)])
	b3_dark_object=min(b3.data[np.where(b3.data>0)])
	b4_dark_object=min(b4.data[np.where(b4.data>0)])
	
	Niwot_utm=[452014,4432155]
	# chimney_utm=[405667,4546640] # old UTM point
	chimney_utm=[405207,4546458]
	if utmpoint==None:
		if niwot:
			UTMpoint=Niwot_utm
		elif chimney:
			UTMpoint=chimney_utm
	radius=5000
		
	xmin=((UTMpoint[0]-radius)-b3.topleft[0])/b3.dx
	ymin=((UTMpoint[1]+radius)-b3.topleft[1])/b3.dy
	xmax=((UTMpoint[0]+radius)-b3.topleft[0])/b3.dx
	ymax=((UTMpoint[1]-radius)-b3.topleft[1])/b3.dy
	
	outputb2=b2.data[ymin:ymax,xmin:xmax].astype('f')
	outputb3=b3.data[ymin:ymax,xmin:xmax].astype('f')
	outputb4=b4.data[ymin:ymax,xmin:xmax].astype('f')
	nodata=np.where((outputb2==0) | (outputb3==0) | (outputb4==0))
	
	outputb2=outputb2*gains[1]+offsets[1]
	outputb3=outputb3*gains[2]+offsets[2]
	outputb4=outputb4*gains[3]+offsets[3]
	
	outputb2=(outputb2-(b2_dark_object*gains[1]+offsets[1]))/Esun[1]
	outputb3=(outputb3-(b3_dark_object*gains[2]+offsets[2]))/Esun[2]
	outputb4=(outputb4-(b4_dark_object*gains[3]+offsets[3]))/Esun[3]
	outputb2[nodata]=0
	outputb3[nodata]=0
	outputb4[nodata]=0
	
	print("Dark Objects:  "+str(b2_dark_object)+'  '+str(b3_dark_object)+'  '+str(b4_dark_object)+"\n")
	
	ndvi=(outputb4-outputb3)/(outputb4+outputb3)
	ndred=(outputb3-outputb2)/(outputb3+outputb2)
	
	res=b2.dx
	origin=np.array(UTMpoint)-np.array([outputb2.shape[1]/2*res,-1*outputb2.shape[0]/2*res])
	return Bunch(ndvi=ndvi,red=ndred,b2=outputb2,b3=b3.data[ymin:ymax,xmin:xmax],b4=outputb4,res=res,origin=origin)
	