#!/usr/bin/env python
import numpy as np
import glob

all_files="wget_complete"
base_filename="wget_example"
var_list=["hus","ta","va","ua"]#,"ps"] #WARNING: not clear how to get ps data, in most cases these all downloaded, or can be read from another file?

# model="ccsm"
# model="cnrm"
# model="miroc"
# model="mri_cgcm3"
model="miroc_esm"

# experiment="historical"
experiment="rcp85"

hist_monthly_req =dict(ccsm=False,cnrm=False,miroc=True,mri_cgcm3=False,miroc_esm=True)
rcp85_monthly_req=dict(ccsm=False,cnrm=False,miroc=True,mri_cgcm3=True, miroc_esm=True)
monthly_req=dict(historical=hist_monthly_req,rcp85=rcp85_monthly_req)[experiment]

hist_urls=dict(ccsm="http://tds.ucar.edu/thredds/fileServer/datazone/cmip5_data/cmip5/output1/"
                +"NCAR/CCSM4/rcp85/6hr/atmos/6hrLev/r6i1p1/v20121031/__VAR__/__NCFILE__",
                cnrm="http://esg.cnrm-game-meteo.fr/thredds/fileServer/esg_dataroot1/CMIP5/output1/"
                +"CNRM-CERFACS/CNRM-CM5/historical/6hr/atmos/6hrLev/r1i1p1/v20110701/__VAR__/__NCFILE__",
                miroc="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MIROC/MIROC5/historical/6hr/atmos/6hrLev/r1i1p1/v20111124/__VAR__/__NCFILE__",
                mri_cgcm3="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MRI/MRI-CGCM3/historical/6hr/atmos/6hrLev/r1i1p1/v20120516/__VAR__/__NCFILE__")

rcp85_urls=dict(mri_cgcm3="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MRI/MRI-CGCM3/rcp85/6hr/atmos/6hrLev/r1i1p1/v20120516/__VAR__/__NCFILE__",
                miroc_esm="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MIROC/MIROC-ESM/rcp85/6hr/atmos/6hrLev/r1i1p1/v20111129/__VAR__/__NCFILE__")
all_urls=dict(historical=hist_urls,rcp85=rcp85_urls)
urls=all_urls[experiment]

Gregorian_object=None #needs to be a class that responds to __getitem__(Year,Month) with ndays... or something like that
calender=dict(ccsm="noleap",cnrm="noleap",miroc="noleap",miroc_esm="gregorian",mri_cgcm3="gregorian")
days_per_month=dict(noleap=[31,28,31,30,31,30,31,31,30,31,30,31],
                   day360=[30]*12,
                   gregorian=Gregorian_object)

dpm=days_per_month[calender[model]]
base_url=urls[model]
generate_monthly=monthly_req[model]
if experiment=="historical":
    start_year=1950
    end_year=2005
else:
    start_year=2006
    end_year=2100
    

def mri_cgcm3_month_based_filename(base,year,month):
    from datetime import datetime, timedelta
    curdate=datetime(year,month,01,00,00)
    if month==12:
        nextdate=datetime(year+1,1,1,00,00)
    else:
        nextdate=datetime(year,month+1,1,00,00)
        
    nextdate-=timedelta(6/24.0) #6hr offset

    greg_date_range="_{0}{1:02}0100-{2}{3:02}{4:02}18.nc\n"
    return base+greg_date_range.format(year,month,nextdate.year,nextdate.month,nextdate.day)
    
def miroc_esm_month_based_filename(base,year,month):
    from datetime import datetime, timedelta
    curdate=datetime(year,month,01,00,00)
    if month==12:
        nextdate=datetime(year+1,1,1,00,00)
    else:
        nextdate=datetime(year,month+1,1,00,00)
        
    greg_date_range="_{0.year}{0.month:02}0106-{1.year}{1.month:02}{1.day:02}00.nc\n"
    return base+greg_date_range.format(curdate,nextdate)


monthly_function=dict(mri_cgcm3=mri_cgcm3_month_based_filename,
                       miroc_esm=miroc_esm_month_based_filename,
                       ccsm=None,miroc=None,cnrm=None)[model]


def create_monthly(filename,template_file):
    """docstring for create_monthly"""
    # strip off the variable name at the beginning e.g. "hus_"
    # and the date range from the end e.g. "_1996010100-1996013118.nc"
    file_base="_".join(template_file.split("_")[1:-1])
    file_base="__VARNAME___"+file_base
    date_range="_{0}{1:02}0100-{0}{1:02}{2}18.nc\n"
    with open(filename,"w") as f:
        for y in range(start_year,end_year+1):
            for m in range(1,13):
                if dpm!=None: # handles most/all noleap and 360day calendar GCMs
                    f.write(file_base+date_range.format(y,m,dpm[m-1]))
                else:
                    f.write(monthly_function(file_base,y,m))

def create_test_file(filename):
    """docstring for create_test_file"""
    max_nfiles=0
    for v in var_list:
        curfiles=glob.glob(v+"_*")
        if len(curfiles)>max_nfiles:
            max_nfiles=len(curfiles)
            curfiles.sort()
            filelist=[f.replace(v,"__VARNAME__") for f in curfiles]
    
    if generate_monthly:
        tmpfile=filelist[0].replace("__VARNAME__",var_list[0])
        create_monthly(filename,tmpfile)
    else:
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
    print("Using {} files out of {}".format(n_out,n_files))
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
    print("""
             -----------------------------------------------------------------------
             Before using the wget-script.sh to get these files, ensure that it has 
             retrieved all the files it can checksum properly by running it until it 
             doesn't get any more.  
             -----------------------------------------------------------------------
             """)
    
    all_files=make_complete_file_list(base_filename)
    unique_filename=make_unique(all_files)
    wget_filename=convert_to_wget_input(unique_filename,base_url)
    
    print("""
             -----------------------------------------------------------------------
             WARNING: ps files are not always generated currently. 
             ps files usually have a different date format which is problematic. 
             Hopefully they were part of the original 1000 files, or can be matched.
             -----------------------------------------------------------------------
             """)
    
    print("")
    print("To download remaining data run : ")
    print("     wget-script.sh -p -F {} &>wget_output_complete &".format(wget_filename))
    print("")
    

if __name__ == '__main__':
    main()