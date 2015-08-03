#!/usr/bin/env python
import sys,os,argparse,traceback
import numpy as np
import glob

all_files="wget_complete"
base_filename="wget_example"
var_list=["hus","ta","va","ua","ps"] #WARNING: not clear how to get ps data, in most cases these all downloaded, or can be read from another file?

hist_urls=dict(ccsm="http://tds.ucar.edu/thredds/fileServer/datazone/cmip5_data/cmip5/output1/"
                +"NCAR/CCSM4/historical/6hr/atmos/6hrLev/r6i1p1/v20121031/__VAR__/__NCFILE__",
                cnrm="http://esg.cnrm-game-meteo.fr/thredds/fileServer/esg_dataroot1/CMIP5/output1/"
                +"CNRM-CERFACS/CNRM-CM5/historical/6hr/atmos/6hrLev/r1i1p1/v20110701/__VAR__/__NCFILE__",
                miroc="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MIROC/MIROC5/historical/6hr/atmos/6hrLev/r1i1p1/v20111124/__VAR__/__NCFILE__",
                mri_cgcm3="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MRI/MRI-CGCM3/historical/6hr/atmos/6hrLev/r1i1p1/v20120516/__VAR__/__NCFILE__",
                miroc_esm="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MIROC/MIROC-ESM/historical/6hr/atmos/6hrLev/r1i1p1/v20111129/__VAR__/__NCFILE__",
                bcc="http://bcccsm.cma.gov.cn/thredds/fileServer/cmip5_data1/output/"
                +"BCC/bcc-csm1-1-m/historical/6hr/atmos/__VAR__/r1i1p1/__NCFILE__")

rcp85_urls=dict(ccsm="http://tds.ucar.edu/thredds/fileServer/datazone/cmip5_data/cmip5/output1/"
                +"NCAR/CCSM4/rcp85/6hr/atmos/6hrLev/r6i1p1/v20121031/__VAR__/__NCFILE__",
                mri_cgcm3="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MRI/MRI-CGCM3/rcp85/6hr/atmos/6hrLev/r1i1p1/v20120516/__VAR__/__NCFILE__",
                miroc_esm="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MIROC/MIROC-ESM/rcp85/6hr/atmos/6hrLev/r1i1p1/v20111129/__VAR__/__NCFILE__",
                miroc="http://dias-esg-nd.tkl.iis.u-tokyo.ac.jp/thredds/fileServer/esg_dataroot/outgoing/output1/"
                +"MIROC/MIROC5/rcp85/6hr/atmos/6hrLev/r1i1p1/v20111124/__VAR__/__NCFILE__",
                cnrm="http://esg.cnrm-game-meteo.fr/thredds/fileServer/esg_dataroot5/CMIP5/output1/"
                +"CNRM-CERFACS/CNRM-CM5/rcp85/6hr/atmos/6hrLev/r1i1p1/v20120525/__VAR__/__NCFILE__",
                hadgem="http://cmip-dn1.badc.rl.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/"
                +"MOHC/HadGEM2-ES/rcp85/6hr/atmos/6hrLev/r1i1p1/v20101208/__VAR__/__NCFILE__")

all_urls = dict(historical=hist_urls,rcp85=rcp85_urls)

hist_monthly_req  = dict(ccsm=False,cnrm=False,miroc=True,mri_cgcm3=False,miroc_esm=False,bcc=False)
rcp85_monthly_req = dict(ccsm=False,cnrm=True,miroc=True,mri_cgcm3=True, miroc_esm=True,bcc=False,hadgem=False)

Gregorian_object=None #Should be a class that responds to __getitem__(Year,Month) with ndays... or something like that
calendar=dict(ccsm="noleap",cnrm="noleap",miroc="noleap",miroc_esm="gregorian",mri_cgcm3="gregorian",bcc="noleap",hadgem="day360")
days_per_month=dict(noleap=[31,28,31,30,31,30,31,31,30,31,30,31],
                   day360=[30]*12,
                   gregorian=Gregorian_object)


def zed_2_18_month_based_filename(base,year,month):
    from datetime import datetime, timedelta
    curdate=datetime(year,month,01,00,00)
    if month==12:
        nextdate=datetime(year+1,1,1,00,00)
    else:
        nextdate=datetime(year,month+1,1,00,00)
        
    nextdate-=timedelta(6/24.0) #6hr offset

    greg_date_range="_{0}{1:02}0100-{2}{3:02}{4:02}18.nc\n"
    return base+greg_date_range.format(year,month,nextdate.year,nextdate.month,nextdate.day)
    
def six_2_zed_month_based_filename(base,year,month):
    from datetime import datetime, timedelta
    curdate=datetime(year,month,01,00,00)
    if month==12:
        nextdate=datetime(year+1,1,1,00,00)
    else:
        nextdate=datetime(year,month+1,1,00,00)
        
    greg_date_range="_{0.year}{0.month:02}0106-{1.year}{1.month:02}{1.day:02}00.nc\n"
    return base+greg_date_range.format(curdate,nextdate)


def global_setup(model="ccsm",experiment="historial"):
    """set up global variables for script
    
    This is not a pretty way to do this, really should return an info struct
    that gets passed around to all routines, but that is a more significant
    refactoring. 
    """
    global monthly_req
    global urls, base_url
    global dpm, generate_monthly
    global start_year, end_year, monthly_function
    
    monthly_req=dict(historical=hist_monthly_req,rcp85=rcp85_monthly_req)[experiment]

    urls=all_urls[experiment]
    base_url=urls[model]

    dpm=days_per_month[calendar[model]]
    if model=="cnrm" and experiment=="rcp85":
        dpm=None
    generate_monthly=monthly_req[model]

    if experiment=="historical":
        start_year=1950
        end_year=2005
    else:
        start_year=2006
        end_year=2100
    
    monthly_function=dict(mri_cgcm3=zed_2_18_month_based_filename,
                           miroc_esm=six_2_zed_month_based_filename,
                           cnrm=six_2_zed_month_based_filename,
                           miroc=None,ccsm=None,bcc=None,hadgem=None)[model]


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

def create_annual(filename,template_file):
    """docstring for create_annual"""
    # strip off the variable name at the beginning e.g. "hus_"
    # and the date range from the end e.g. "_1996010100-1996013118.nc"
    file_base="_".join(template_file.split("_")[1:-1])
    file_base="__VARNAME___"+file_base
    date_range="_{0}010100-{0}123118.nc\n"
    with open(filename,"w") as f:
        for y in range(start_year,end_year+1):
            f.write(file_base+date_range.format(y))

def create_test_file(filename,annual_ps):
    """docstring for create_test_file"""
    max_nfiles=0
    for v in var_list:
        curfiles=glob.glob(v+"_*")
        if len(curfiles)>max_nfiles:
            max_nfiles=len(curfiles)
            curfiles.sort()
            filelist=[f.replace(v,"__VARNAME__") for f in curfiles]
    if annual_ps:
        print("Creating an annual list")
        tmpfile=filelist[0].replace("__VARNAME__",var_list[0])
        create_annual(filename,tmpfile)
    elif generate_monthly:
        print("Creating a monthly list")
        tmpfile=filelist[0].replace("__VARNAME__",var_list[0])
        create_monthly(filename,tmpfile)
    else:
        print("Creating file with {} files per variable.".format(max_nfiles))
        with open(filename,"w") as f:
            for l in filelist:
                f.write(l+"\n")
    
    

def make_complete_file_list(filename,annual_ps):
    """create a file with all dates and variables"""
    complete_file=all_files
    if not glob.glob(filename):
        print("Creating file...")
        create_test_file(filename,annual_ps)
    else:
        print("Using existing file:"+filename)
    with open(filename,"rU") as f:
        fout=open(complete_file,"w")
        for l in f:
            if annual_ps:
                fout.write(l.replace("__VARNAME__","ps"))
            else:
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
            varname=l.split("_")[0] # filenames are e.g. "hus_6hr_..."
            ncfilename=l.strip()    # remove trailing newline
            cur_url=url.replace("__VAR__",varname).replace("__NCFILE__",ncfilename)
            fout.write(outputline.format(ncfilename,cur_url))
        fout.close()
    return wget_file

def make_fill_file(filename):
    """Make a file with a list of files present in the current directory with 0 size
    
    equivalent to : ls -lh *.nc | grep " 0 " | awk '{ print $9}' >filename+"_fill"
    """
    files=glob.glob("*.nc")                     # get a list of relevant files
    with open(filename+"_fill","w") as fout:    # open the output file
        for f in files:                         # loop through files
            fileinfo = os.stat(f)               # get file size information
            if fileinfo.st_size==0:             # if file size is 0 then
                fout.write(f+"\n")                 # write the filename to the output file
    return filename+"_fill"
    

def main(model="ccsm",experiment="historical",base=None,all_files=None, unique_filename=None, fill_missing=False,annual_ps=False):
    """Set up an input file for wget"""
    print("""
             -----------------------------------------------------------------------
             Before using the wget-script.sh to get these files, ensure that it has 
             retrieved all the files it can checksum properly by running it until it 
             doesn't get any more.  
             -----------------------------------------------------------------------
             """)
    
    try:
        global_setup(model=model,experiment=experiment)
    except KeyError:
        traceback.print_exc()
        sys.exit("script not set up to work for model {0}, experiment {1}".format(model,experiment))
        
    
    if (unique_filename==None) and (not fill_missing):
        if all_files==None:
            all_files=make_complete_file_list(base_filename,annual_ps)
        unique_filename=make_unique(all_files)
    if (fill_missing==True):
        print("Filling")
        unique_filename=make_fill_file(base_filename)
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
    try:
        parser= argparse.ArgumentParser(description='Generate a wget inputfile to additional CMIP5 data. ')
        parser.add_argument('model',nargs="?",default="ccsm",action='store',
                            help="Model to process [ccsm,cnrm,mri_cgcm3,miroc,miroc_esm,etc]")
        parser.add_argument('experiment',nargs="?",action='store', default="historical",
                            help="Experiment to process <historical,rcp85>")
        parser.add_argument('-b','--base',dest='base_filename',nargs="?",action='store',default=None,
                            help="Base filename to use in creating other files")
        parser.add_argument('-a','--all',dest='all_files',nargs="?",action='store',default=None,
                            help="Filename of file containing all files to test for downloading (if a file exists it will be skipped).")
        parser.add_argument('-u','--unique',dest='unique_file',nargs="?",action='store',default=None,
                            help="Filename of file containing list of files to download")
        parser.add_argument('-v', '--version',action='version',
                            version='get_remaining_files v1.0')
        parser.add_argument('-f','--fill',action='store_true',
                            default=False, dest="fill",help="Attempt to redownload zero length files in current directory.")
        parser.add_argument('-p','--print',action='store_true',
                            default=False, dest="print_models",help="Print available internal model names.")
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        parser.add_argument ('--annual_ps', action='store_true',
                default=False, help='write ps files as annual', dest='annual_ps')
        args = parser.parse_args()
        
        if args.print_models:
            print(calendar.keys())
            sys.exit(0)
        
        exit_code = main(model=args.model,
                         experiment=args.experiment,
                         base=args.base_filename,
                         all_files=args.all_files,
                         unique_filename=args.unique_file, 
                         fill_missing=args.fill,
                         annual_ps=args.annual_ps)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
