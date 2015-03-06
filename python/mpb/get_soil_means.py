#!/usr/bin/env python

from load_data import cols as load_cols
import date_fun
import numpy as np
import iterdim
import glob

def load_data(filename):
	
	f=open(filename,'r')
	for i in range(3):l=f.next()
	d=np.zeros((100000,200))
	maxwidth=0
	i=0
	for l in f:
		curdata=np.array(l.split(),'d')
		maxwidth=np.max([maxwidth,curdata.size])
		d[i,0:curdata.size]=curdata
		i+=1
	return d[0:i,0:maxwidth]


def get_stats(data,minmax=[0,0],qccheck=None):
	output=np.zeros((len(data[:,0]),3))
	print(minmax)
	for i in range(output.shape[0]):
		if qccheck==None:
			tmp=np.where((data[i,:]>minmax[0])&(data[i,:]<minmax[1]))
		else:
			tmp=np.where((data[i,:]>minmax[0])&(data[i,:]<minmax[1])&(qccheck[i,:]==1))
		if len(tmp[0])>0:
			output[i,0]=np.mean(data[i,:][tmp])
			output[i,1]=np.min(data[i,:][tmp])
			output[i,2]=np.max(data[i,:][tmp])
		else:
			output[i,:]=-9999
	return output


def write_file(name,dates,data):
	output=np.hstack((dates,data))
	np.savetxt(name+'_stats.txt',output,fmt=('%i '*6)+('%f '*3))


def compile_files(filesearch,dt=1.0/24/4.0):
	# find all files to compile
	files=glob.glob(filesearch)
	# make list variables to stores data and dates
	alld=list()
	mjds=list()
	# set up dummy min and max dates to store real min and max dates
	minmjd=date_fun.date2mjd(2050,1,1,1,1,1)
	maxmjd=date_fun.date2mjd(1900,1,1,1,1,1)
	# start with 6 columns (for the date) and add more as we go
	ncols=6
	# loop over files reading in data
	for f in files:
		alld.append(load_cols(f))
	# now loop over data getting dates and column numbers in each
	for d in alld:
		mjds.append(np.round(date_fun.datearr2mjd(d[:,0:6])*4*24)/4.0/24.0)
		# number of columns in this file not including the date columns
		ncols+=d[0,6:].size
	# loop over all date arrays looking for min and max dates 
	for mjd in mjds:
		minmjd=np.min((minmjd,np.min(mjd)))
		maxmjd=np.max((maxmjd,np.max(mjd)))
		
	# assumes the time step is 15minutes
	# dt=1.0/24/4.0
	# calculate the number of time steps in the output data
	ntimes=(maxmjd-minmjd)/dt+2
	# create an array of dates
	outputmjd=np.arange(minmjd,maxmjd+dt,dt)
	# create the output data array
	outputdata=np.zeros((ntimes,ncols))
	# add the dates to the output data array
	outputdata[:,0:6]=date_fun.mjd2date(outputmjd,roundseconds=True)
	
	# now the "hard" part, loop over all data and 
	# then put that data in its place in the output data array
	startcol=6
	for (d,mjd) in zip(alld,mjds):
		curcols=d[0,6:].size
		for (curd,curm) in zip(iterdim.iterdim(d),mjd):
			curloc=np.round((curm-minmjd)/dt)
			outputdata[curloc,startcol:startcol+curcols]=curd[6:]
		startcol+=curcols
	np.savetxt('combined.txt',outputdata,fmt=('%i '*6)+'%f '*(ncols-6))
	

def get_soil_means(filename):
	data=load_data(filename)
	indices=[10,30,60,-10,-30,-60]
	minmax=np.array([[0.05,0.45],[-30,60]])
	names=['smc10','smc30','smc60','ts10','ts30','ts60']
	dates=data[1:,0:6]
	# if (data[0,0]<39000):
	# 	dates=date_fun.mjd2date(date_fun.macexcel2mjd(data[1:,0]),roundseconds=True)
	# else:
	# 	dates=date_fun.mjd2date(date_fun.excel2mjd(data[1:,0]),roundseconds=True)
	
	for i in range(len(indices)):
		columns=np.where(data[0,:]==indices[i])[0]
		qc=columns+1
		curdata=get_stats(data[1:,columns],minmax=minmax[indices[i]<0,:])
		# if indices[i]<0:
		# 	curdata=get_stats(data[1:,columns],minmax=minmax[indices[i]<0,:])
		# else:
		# 	curdata=get_stats(data[1:,columns],minmax=minmax[indices[i]<0,:],qccheck=data[1:,qc])
		write_file(names[i],dates,curdata)

