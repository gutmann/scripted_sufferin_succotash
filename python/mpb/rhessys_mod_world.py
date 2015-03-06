#!/usr/bin/env python
# encoding: utf-8
"""
rhessys_mod_world.py

Created by Ethan Gutmann on 2011-06-28.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""

import sys
import getopt
import re

import numpy as np
import load_data


help_message = '''
rhessys_mod_world

USAGE
	rhessys_mod_world <worldfile> <variable> <mapfile> <clumpfile> 

Modifies a RHESSys worldfile to update "variable" with values stored in mapfile.
Assumes clumpfile and mapfile are text file grids that overlap each other perfectly. 

The current worldfile is saved as <worldfile>.pre_<variable>_update

EXAMPLE:
	rhessys_mod_world niwot_2d_worldfile cs_leafc newLAI.txt niwot.cl
	
BUGS:
	I'm not sure if this is a bug or not (I don't know how rhessys works internally)
	If there are two patches with the same patch ID, but that are in different zones
	It assigns the mean for the combined patches to both patches.  It does not
	split patches based on their zone. This should not have much effect, but it is 
	present. 
'''


class Usage(Exception):
	
	def __init__(self, msg):
		self.msg = msg
	


def update_worldfile(worldfile,variable, mapfile,clumpfile,gain=False,delta=False):
	mapdata=load_data.cols(mapfile)
	clumps=load_data.cols(clumpfile)
	
	# set the inital state to be out of a patch
	inpatch=False
	
	# open the outputfile
	fo=open('tempworldfile_'+variable,'w')
	#loop through the input file replacing values as feasible
	with open(worldfile,'r') as f:
		for line in f: 
			# if we have already found that we are within a valid clump/patch
			if inpatch:
				# search for the variable of interest, if it matches, replace
 				# the value and set the state as out of the clump/patch
				if re.match('^ *-*[0-9]?\.*[0-9]* *'+variable+' *$'):
					if delta:
						thisvalue=np.float(line.split()[0])
						fo.write('                '+str(curvalue+thisvalue)+'     '+variable+'\n')
					elif gain:
						thisvalue=np.float(line.split()[0])
						fo.write('                '+str(curvalue*thisvalue)+'     '+variable+'\n')
					else:
						fo.write('                '+str(curvalue)+'     '+variable+'\n')
						
					inpatch=False
				else:
					# if we are not at the variable of interest just write the
					# current line to the outputfile
					fo.write(line)
			# if we are on a line that looks like a patch_ID:
			elif re.match('^ *[0-9]? *patch_ID *$'):
				# get the patchID
				patchID=np.int(line.split()[0])
				# find the patchID locations in the clump file
				tmp=where(clumps==patchID)
				# calculate the mean value at those locations in the mapfile
				curvalue=np.mean(mapdata[tmp])
				# and write the current line to the output file
				fo.write(line)
				# then set the state to be in a patch
				inpatch=True
			else:
				# else we didn't match any special cases, just write the 
				# current line to the output file
				fo.write(line)
			
		
	# close the output file (input file is automagically closed)
	fo.close()


# Just the main program, parse the commandline arguments and pass the to
# update_worldfile (the main worker procedure)
def main(argv=None):
	delta=None
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
			if option in ("--delta","-d")
				delta=True
			if option in ("--gain","-g")
				gain=True
		if len(args) == 4:
			update_worldfile(args[0],args[1],args[2],args[3],delta=delta,gain=gain)
		else: 
			raise Usage(msg)
	
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2
		
	return 0


if __name__ == "__main__":
	sys.exit(main())
