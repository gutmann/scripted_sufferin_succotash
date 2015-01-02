#!/usr/bin/env python
# encoding: utf-8
"""
rhessys_changeXinY.py

Created by Ethan Gutmann on 2011-06-10.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""

import sys
import getopt
import re

help_message = '''
rhessys_changeXinY [-mul] [-add] <filename> <variable> <value> [searchfield] [newfilename]

Reads through <filename> searching for <variable> (e.g. cs_leafc) and replacing
it's current value with <value>.  

If searchfield is given, it first searches for this field, after finding the 
searchfield it replaces the first instance of <variable> then searches for the
searchfield again.	Searchfield can include matches for multiple lines, these
lines should be split by the two characters \$

searchfield and variable can both be regular expressions

If newfilename is specified, it stores the results in a file with that name, 
otherwise it adds "variable=value" to the original filename.  

if -mul set, multiplies the existing value by the new value
if -add set, adds the existing value to the new value

'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def findfield(fin,field,fout):
	'''search fin file handle for the regular expression field
	write all lines that donot match to fout and return the line that matches.
	'''
	line=fin.next()
	while (line != None) and (not re.match(field,line)):
		fout.write(line)
		try:
			line=fin.next()
		except StopIteration:
			line=None
	return line

def replace(filename,variable,value,searchfield=None,outputfile=None,mul=None,add=None):
	
	if outputfile==None:
		outputfile=filename+'_'+str(variable)+'_'+str(value)
	
	if searchfield!=None:
		searchfield=searchfield.split('\$')
		print(repr(searchfield))
	
	# open the input and output files
	fout=open(outputfile,'w')
	with open(filename,'r') as f:
		curline=''
		while curline !=None:
			# if a search field was specified find it (or them) first
			if searchfield!=None:
				# loop over all searchfields
				for searchable in searchfield:
					# if curline is None then we hit the end of the file already
					if curline != None:
						# find the current searchfield anywhere on a line
						curline=findfield(f,'.*'+searchable+'.*$',fout)
						# if we found it, print it to the new file
						if curline !=None:
							fout.write(curline)
							
			# now find the variable we are replacing
			if curline != None:
				# find the current variable name, it must be followed by nothing 
				# but whitespace on the current line and have at least one space
				# in front of it. 
				curline=findfield(f,'.* +'+variable+' *$',fout)
			if curline != None:
				# and write it to the new file with the new value
				newline=curline.split()
				altsplit=curline.split(newline[0])
				if mul!=None:
					curvalue=float(newline[0])*float(value)
				elif add!=None:
					curvalue=float(newline[0])+float(value)
				else:
					curvalue=value
				fout.write(altsplit[0]+str(curvalue)+altsplit[1])
				
	fout.close


def main(argv=None):
	'''The main program just handles command line options and passes them to replace''' 
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "maho:v", ["mul","add","help", "output="])
		except getopt.error, msg:
			raise Usage(msg)
		mul=None
		add=None
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-m", "--mul"):
				mul=1
			if option in ("-a", "--add"):
				add=1
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-o", "--output"):
				output = value
	
		if len(args)<3:
			raise Usage(help_message)
		elif len(args)==3:
			replace(args[0],args[1],args[2],mul=mul,add=add)
		elif len(args)==4:
			replace(args[0],args[1],args[2],searchfield=args[3],mul=mul,add=add)
		elif len(args)==5:
			replace(args[0],args[1],args[2],searchfield=args[3],outputfile=args[4],mul=mul,add=add)
			
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	


if __name__ == "__main__":
	sys.exit(main())
