#!/usr/bin/env python

import glob
import mygis

vars2load=[""]

def load_data(filename):
    

def convert2daily(filename):
    print("  Loading")
    forcing = load_data(filename)
    
    print("  Converting")
    daily = convert_data(forcing)
    
    print("  Writing")
    write_daily(daily, filename)
    

def main():
    files=glob.glob("icar_*_00")
    files.sort()
    
    for f in files:
        print(f)
        convert2daily(f)

if __name__ == '__main__':
    main()
