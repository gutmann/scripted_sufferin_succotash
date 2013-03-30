'''
Created on Mar 31, 2011

@author: gutmann

	date_fun
		a collection of functions for working with oddball (I'm looking at you 
		Excel) dates.  
		
		In particular
		excel2julian - convert 1899-12-30 based dates to julian dates
		macexcel2julian - convert 1904-1-1 based dates to julian dates
			2date - convert to a list (year, month,day,hour,minute,second)
			2mjd  - convert to modified Julian Day
		
'''
import julday
from numpy import double,zeros,floor,ones,array
import numpy as np

def excel2julian(date_xls):
	if np.array(date_xls).size==1:
		return julday.jul_day(1899,12,30,0,0,0)+date_xls
	
	output=np.zeros(np.array(date_xls).size)
	for i in range(np.array(date_xls).size):
		output[i]=julday.jul_day(1899,12,30,0,0,0)+date_xls[i]
	return output

def macexcel2julian(date_xls):
	if np.array(date_xls).size==1:
		return julday.jul_day(1904,1,1,0,0,0)+date_xls-1
	
	output=np.zeros(np.array(date_xls).size)
	for i in range(np.array(date_xls).size):
		output[i]=julday.jul_day(1904,1,1,0,0,0)+date_xls[i]-1
	return output

def excel2mjd(date_xls):
	if np.array(date_xls).size==1:
		return julday.mjul_day(1899,12,30,0,0,0)+date_xls
	
	output=np.zeros(np.array(date_xls).size)
	for i in range(np.array(date_xls).size):
		output[i]=julday.mjul_day(1899,12,30,0,0,0)+date_xls[i]
	return output

def macexcel2mjd(date_xls):
	if np.array(date_xls).size==1:
		return julday.mjul_day(1904,1,1,0,0,0)+date_xls-1
	
	output=np.zeros(np.array(date_xls).size)
	for i in range(np.array(date_xls).size):
		output[i]=julday.mjul_day(1904,1,1,0,0,0)+date_xls[i]-1
	return output

def excel2date(date_xls):
	if np.array(date_xls).size==1:
		return julday.caldate(julday.jul_day(1899,12,30,0,0,0)+date_xls)
	
	dates=[julday.caldate(julday.jul_day(1899,12,30,0,0,0)+date_xls[i]) for i in range(np.array(date_xls).size)]
	return dates
	

def macexcel2date(date_xls):
	if np.array(date_xls).size==1:
		return julday.caldate(julday.jul_day(1904,1,1,0,0,0)+date_xls-1)
	
	dates=[julday.caldate(julday.jul_day(1904,1,1,0,0,0)+date_xls[i]-1) for i in range(np.array(date_xls).size)]
	return dates
	

def mjd2datetime(mjd,roundseconds=None):
	return julday2datetime(double(mjd)+julday.MJD0,roundseconds=roundseconds)

def mjd2date(mjd,roundseconds=None):
	return julday2date(double(mjd)+julday.MJD0,roundseconds=roundseconds)


def date2mjd(year,month,day,hour,minute,second=None):
	year=array(year)
	ndates=year.size
	
	if second==None: second=zeros(ndates)
	mjd=zeros(ndates)
	if ndates==1:return julday.mjul_day(np.float(year),month,day,hour,minute,np.float(second))
	for i in np.arange(ndates):
		mjd[i]=julday.mjul_day(year[i],month[i],day[i],hour[i],minute[i],second[i])
	return mjd 
	
def date2julian(year,month,day,hour,minute,second=None):
	try:
		ndates=len(year)
	except TypeError:
		ndates=1
	
	if second==None: second=zeros(ndates)
	jd=zeros(ndates)
	for i in range(ndates):
		jd[i]=julday.jul_day(year[i],month[i],day[i],hour[i],minute[i],second[i])
	return jd

def ydoyhm2mjd(year,doy,hhmm):
	return ydoyhm2julday(year,doy,hhmm)-julday.MJD0

def ydoyhm2julday(inyear,doy,hhmm):
	year=array(inyear)
	hours=floor(hhmm/100)
	minutes=(hhmm%100)
	if year.size > 1:
		jdays=array([julday.jul_day(year[i],1,1,hours[i],minutes[i],0) for i in range(year.size)])
	else:
		jdays=julday.jul_day(inyear,1,1,hours,minutes,0)
	return jdays+doy

def julday2datetime(jdayin,roundseconds=None):
	import datetime
	jday=np.array(jdayin)
	if jday.size==1:
		if roundseconds==None:
			(year,month,day,hour,minute,second)=julday.caldate(jdayin)
		else:
			(year,month,day,hour,minute,second)=julday.caldate(jdayin+0.00034)
			second=0
		return datetime.datetime(int(year),int(month),int(day),int(hour),int(minute),int(round(second))%60)
	
	datetimes=[datetime.datetime(1990,1,1,0,0,0) for i in range(jday.size)]
	for i in range(len(jday)):

		if roundseconds == None:
			(year,month,day,hour,minute,second)=julday.caldate(jday[i])
		else:
			(year,month,day,hour,minute,second)=julday.caldate(jday[i]+0.00034)
			second=0

#		 import pdb; pdb.set_trace()
		datetimes[i]=datetime.datetime(int(year),int(month),int(day),int(hour),int(minute),int(round(second))%60)
		
	return datetimes

def julday2date(jdayin,roundseconds=None):
	jday=np.array(jdayin)
	if jday.size==1:
		if roundseconds==None:
			return julday.caldate(jdayin)
		else:
			(year,month,day,hour,minute,second)=julday.caldate(jday+0.00034)
			second=0
			return np.array([year,month,day,hour,minute,second])
	
	output=np.zeros((jday.size,6))
	for i in range(len(jday)):
		
		if roundseconds == None:
			(year,month,day,hour,minute,second)=julday.caldate(jday[i])
		else:
			(year,month,day,hour,minute,second)=julday.caldate(jday[i]+0.00034)
			second=0
		output[i,:]=np.array([year,month,day,hour,minute,second])
	return output


def datetime2julday(dt):
	juldays=np.zeros(len(dt),dtype='d')
	for i in range(len(juldays)):
		juldays[i]=julday.jul_day(dt[i].year,dt[i].month,dt[i].day,dt[i].hour,dt[i].minute,dt[i].second)
		
	return juldays

def datetime2mjd(dt):
	return datetime2julday(dt)-julday.MJD0


def datearr2mjd(dates):
	
	sz=dates.shape
	if len(sz)==1:
		if sz[0]==3:
			return date2mjd(dates[0],dates[1],dates[2],12,0)
		if sz[0]==4:
			return date2mjd(dates[0],dates[1],dates[2],dates[3],0)
		if sz[0]==5:
			return date2mjd(dates[0],dates[1],dates[2],dates[3],dates[4])
		if sz[0]==6:
			return date2mjd(dates[0],dates[1],dates[2],dates[3],dates[4],second=dates[5])
	else:
		if sz[1]==3:
			return date2mjd(dates[:,0],dates[:,1],dates[:,2],np.zeros(sz[0])+12,np.zeros(sz[0]))
		if sz[1]==4:
			return date2mjd(dates[:,0],dates[:,1],dates[:,2],dates[:,3],np.zeros(sz[0]))
		if sz[1]==5:
			return date2mjd(dates[:,0],dates[:,1],dates[:,2],dates[:,3],dates[:,4])
		if sz[1]==6:
			return date2mjd(dates[:,0],dates[:,1],dates[:,2],dates[:,3],dates[:,4],second=dates[:,5])

	
def datearr2datetime(dates):
	import datetime
	sz=dates.shape
	if sz[1]==0:
		if sz[0]==3:
			hours=12
		else:
			hours=dates[3]
		if sz[0]<=4:
			minutes=0
		else:
			minutes=dates[4]
		if sz[0]<=5:
			seconds=0
		else:
			seconds=dates[5]
		return datetime.datetime(dates[0],dates[1],dates[2],hours,minutes,seconds)
	else:
		if sz[1]==3:
			hours=np.zeros(sz[0])+12
		else:
			hours=dates[:,3]
		if sz[1]<=4:
			minutes=np.zeros(sz[0])
		else:
			minutes=dates[:,4]
		if sz[1]<=5:
			seconds=np.zeros(sz[0])
		else:
			seconds=dates[:,5]

		tmp=np.where(dates[:,0]<1950)
		if tmp[0].size>0:
			dates[tmp,0]=1950
			dates[tmp,1]=1
			dates[tmp,2]=1
		datetimes=[datetime.datetime(dates[i,0].astype('i'),dates[i,1].astype('i'),dates[i,2].astype('i'),hours[i].astype('i'),minutes[i].astype('i'),seconds[i].astype('i')) for i in range(sz[0])]
		
	return datetimes
