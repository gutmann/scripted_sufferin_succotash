#!/usr/bin/env python

import glob,os,re
from math import floor


def get_year(filename):
    """docstring for get_year"""
    with open(filename,"Ur") as f:
        for l in f:
            if re.match(".*output_file.*",l):
                year=int(l.split("_")[-1][:4])
    
    return year

def main():
    """docstring for main"""
    options_file=glob.glob("icar_*opt.nml")[0]
    cur_year=get_year(options_file)
    
    files=glob.glob("output/icar_"+str(cur_year)+"-0*")
    files.sort()
    last=files[-1]
    last_number=int(last.split("-")[-1])
    if last_number>=8760:
        # last_number=8760
        restart_step=1
        restart_file="output/icar_"+str(cur_year)+"-08760"
        cur_year+=1
    else:
        restart_step=int(floor(last_number/6)+1)
        restart_file="output/icar_"+str(cur_year)+"-{:05}".format(restart_step*6-6)
    
    # escape the filename for sed
    restart_file=restart_file.replace("/","\/")
    os.system("sed -e 's/__RESTART_FILE__/"+restart_file   +"/g;"+
                      "s/__RESTART_STEP__/"+str(restart_step)   +"/g;"+
                      "s/__YEAR__/"        +str(cur_year)  +"/g;"+
                      "s/__NEXTYEAR__/"    +str(cur_year+1)+"/g;"+
                      "' template.nml >"+options_file)
    
        
    
    

if __name__ == '__main__':
    main()