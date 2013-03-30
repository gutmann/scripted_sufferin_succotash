#!/Library/Frameworks/Python.framework/Versions/2.6/bin/python
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
		
# import sys
# print is_number(sys.argv[1])
