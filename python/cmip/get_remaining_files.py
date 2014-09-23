#!/usr/bin/env python
import numpy as np
import glob

all_files="wget_complete"
base_filename="wget_example"
var_list=["hus","ta","va","ua"]#,"ps"] #WARNING: not clear how to get ps data, in most cases these all downloaded, or can be read from another file?

# CCSM data
base_url="http://tds.ucar.edu/thredds/fileServer/datazone/cmip5_data/cmip5/output1/NCAR/CCSM4/rcp85/6hr/atmos/6hrLev/r6i1p1/v20121031/__VAR__/__NCFILE__"

def create_test_file(filename):
    """docstring for create_test_file"""
    max_nfiles=0
    for v in var_list:
        curfiles=glob.glob(v+"_*")
        if len(curfiles)>max_nfiles:
            max_nfiles=len(curfiles)
            curfiles.sort()
            filelist=[f.replace(v,"__VARNAME__") for f in curfiles]
    
    print("Creating file with {} files per variable.".format(max_nfiles))
    with open(filename,"w") as f:
        for l in filelist:
            f.write(l+"\n")
    
    

def make_complete_file_list(filename):
    """create a file with all dates and variables"""
    complete_file=all_files
    if not glob.glob(filename):
        create_test_file(filename)
    with open(filename,"rU") as f:
        fout=open(complete_file,"w")
        for l in f:
            for v in var_list:
                fout.write(l.replace("__VARNAME__",v))
        fout.close()
        
    return complete_file

def make_unique(filename):
    """Select files from input file data that don't exist in the current directory"""
    n_files=0
    n_out=0
    unique_filename=filename+"_unique"
    with open(filename,"rU") as f:
        fout=open(unique_filename,"w")
        for l in f:
            n_files+=1
            if glob.glob(l.strip()):
                pass
            else:
                n_out+=1
                fout.write(l)
        fout.close()
    print("Found {} files out of {}".format(n_out,n_files))
    return unique_filename

def convert_to_wget_input(filename,url):
    """convert a list of files to be downloaded into a wget-script inputfile"""
    wget_file=filename+"_wget"
    outputline="'{0}' '{1}' 'MD5' '00000000000000000000000000000000'\n"
    with open(filename,"rU") as f:
        fout=open(wget_file,"w")
        for l in f:
            varname=l.split("_")[0] #filenames are e.g. "hus_6hr_..."
            ncfilename=l.strip() #remove trailing newline
            cur_url=url.replace("__VAR__",varname).replace("__NCFILE__",ncfilename)
            fout.write(outputline.format(ncfilename,cur_url))
        fout.close()
    return wget_file

def main():
    """Set up an input file for wget"""
    print("""Before using the wget-script.sh to get these files, ensure that it has 
             retrieved all the files it can checksum properly by running it until it 
             doesn't get any more.  """)
    
    all_files=make_complete_file_list(base_filename)
    unique_filename=make_unique(all_files)
    wget_filename=convert_to_wget_input(unique_filename,base_url)
    
    print("""WARNING: ps files are not generated currently. 
        ps files usually have a different date format so this script doesn't work. 
        Hopefully they were part of the original 1000 files...""")
    
    print("")
    print("To download remaining data run : ")
    print("     wget-script.sh {} &>wget_output_complete &".format(wget_filename))
    

if __name__ == '__main__':
    main()