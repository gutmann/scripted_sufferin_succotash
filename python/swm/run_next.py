#!/usr/bin/env python
import glob
import os
import sys
import re
import subprocess

def get_last_file(output_dir):
    """docstring for get_last_file"""
    files=glob.glob(output_dir+"swim_out[0-9]*")
    files.sort()
    return files[-1]
    
def get_next_time(full_path):
    """docstring for fname"""
    base_filename="swim_out"
    filename=full_path.split("/")[-1]
    next_time=int(filename[len(base_filename):])+1
    return next_time

def setup_file(last_file,next_time):
    """docstring for setup_file"""
    filename="real_options.namelist"
    with open(filename,"ru") as f:
        fout=open(filename+".temp","w")
        for l in f:
            if re.match("\s*restart\s*=.*",l):
                l="restart=true,\n"
            if re.match("\s*restart_step\s*=.*",l):
                l="restart_step={},\n".format(next_time)
            if re.match("\s*restart_file\s*=.*",l):
                l="restart_file=\"{}\"\n".format(last_file)
            
            fout.write(l)
        fout.close()
    
    os.rename(filename,filename+".last")
    os.rename(filename+".temp",filename)

def launch():
    retcode=subprocess.call("bsub",stdin=open("batch_submit.sh"))

def main(output_dir="output/"):
    """docstring for main"""
    
    last_file=get_last_file(output_dir)
    next_time=get_next_time(last_file)
    
    print(last_file)
    print(next_time)
    setup_file(last_file,next_time)
    launch()
    
def usage():
    """docstring for usage"""
    print("""Usage: run_next.py [-h|--help] [output_dir]
    
             output_dir default = 'output_dir/'
             -h|--help print this message""")

if __name__ == '__main__':
    args=sys.argv
    if len(args)>1:
        if (args[1]=="-h") or (args[1]=="--help"):
            usage()
            os._exit(0)
        else:
            output_dir=args[1]
    else:
        output_dir="output/"
    
    main(output_dir=output_dir)