#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Ethan Gutmann on 2011-12-01.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""

import sys
import getopt


help_message = '''
convert RHESSys daily output files to monthly averages

daily2monthly <dailyfilename> <monthlyfilename>
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

import load_data
import numpy as np

def daily2monthly(dailyfilename, monthlyfilename):
	d=load_data.cols(dailyfilename)
	outputdata=d.copy()
	i=0
	iout=0
	nout=0
	while (i<len(d[:,0])):
		curyear=d[i,2]
		curmonth=d[i,1]
		outputdata[iout,0:3]=np.array([0,curmonth,curyear])
		nout=0
		in_month=True
		while (i<len(d[:,0])) & in_month:
			in_month=((d[i,2] == curyear) & (d[i,1]==curmonth))
			if in_month:
				outputdata[iout,3:]+=d[i,3:]
				nout+=1
				i+=1
		# convert accumulations to averages
		outputdata[iout,3:]/=nout
		# convert water balance terms back to accumulations instead of averages. 
		outputdata[iout,[5,10,11,12,13,15,16,17,18,30,34]]*=nout
		iout+=1
	
	finaloutput=outputdata[0:iout,:]
	np.savetxt(monthlyfilename,finaloutput,fmt='%f')
	

def main(argv=None):
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-o", "--output"):
				output = value
		if len(args)<2:
			raise Usage(help_message)
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2
	daily2monthly(args[0],args[1])


if __name__ == "__main__":
	sys.exit(main())
