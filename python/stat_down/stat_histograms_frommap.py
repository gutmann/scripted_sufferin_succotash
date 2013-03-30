#!/usr/bin/env python

"""
SYNOPSIS

    template_argparse.py [-h] [--verbose] [-v, --version] <filename>

DESCRIPTION

    TODO This describes how to use this script.
    This docstring will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION

    
"""

import sys
import os
import traceback
import argparse
import glob
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import swim_io as mygis

def plot_for(filename=None,files=None,varname="data",plotzero=False, trailing="n",nbins=50,data=None,
             xlim=None,xlabel=None,title=None,figname="current_histogram.png",ylim=None):

    if files==None:
        files=[f for f in glob.glob(filename) if (not re.match(".*dry_"+trailing+".*",f) 
                                            and not re.match(".*wet_"+trailing+".*",f) 
                                            and not re.match(".*hot_"+trailing+".*",f) 
                                            and not re.match(".*cold_"+trailing+".*",f))]
    if data==None:
        data=[]
        vmin=9999
        vmax=-9999
        for f in files:
            data.append(mygis.read_nc(f,varname).data)
            vmin=min(vmin,data[-1].min())
            vmax=max(vmax,data[-1].max())
        
    hdata=[]
    if xlim==None:
        xlim=[vmin,vmax]
    xlim=np.array(xlim,dtype="f")
    if re.match(".*_pr_RMS_map.nc",files[0]): xlim=xlim/365
    for d in data:
        if re.match(".*_pr_RMS_map.nc",files[0]): d/=365.0
        hdata.append(np.histogram(d,bins=nbins,range=xlim))
    plt.clf()
    for h,f in zip(hdata,files):
        x=np.array((h[1][:-1]+h[1][1:]))/2.0
        y=h[0]
        color="red"
        linewidth=2.0
        linestyle="-"
        marker=None
        if re.match("CA.*",f) or re.match(".*_CA.*",f): color="blue"
        if re.match("ncep_.*",f) or re.match(".*_ncep_.*",f): color="dark"+color
        if re.match("SDmon_.*",f) or re.match(".*_SDmon_.*",f): linestyle="--"
        if re.match(".*e0_"+trailing+".*",f): linestyle=":"
        if re.match(".*SARe0_"+trailing+".*",f): 
            linestyle="-"
            color="green"
        if re.match(".*e1_"+trailing+".*",f): linestyle="--"
        if re.match("obs.*",f): 
            color="black"
            linewidth=3.0
            marker="o"
        # if re.match(".*_pr_RMS_map.nc",f):x/=365.0
        
        if (not re.match(".*e_"+trailing+".*",f) 
            and not re.match(".*Ae0_"+trailing+".*",f)
            and not re.match(".*De0_"+trailing+".*",f)
            and not re.match(".*e1_"+trailing+".*",f)):
            plt.plot(x,y,color=color,marker=marker,linewidth=linewidth,linestyle=linestyle,
                    label=" ".join(f.split(".")[0].split("_")[:2]))
        
    plt.ylabel("N")
    if ylim:plt.ylim(ylim)
    if xlabel:plt.xlabel(xlabel)
    if plotzero:
        plt.plot([0,0],plt.gca().get_ylim(),color="black",linewidth=1)
    plt.legend(loc=0,ncol=2,prop=matplotlib.font_manager.FontProperties(size="small"))
    # plt.ylim(0,20000)
    if xlim!=None:plt.xlim(xlim)
    print(figname)
    plt.savefig(figname)
    
def stat_v_obs():  
    filesearches=["*tasmax_bias_map.nc",
                  "*tasmax_RMS_map.nc",
                  "*tasmin_bias_map.nc",
                  "*tasmin_RMS_map.nc",
                  "*pr_bias_map.nc",
                  "*pr_RMS_map.nc"]
    plotzeros=[True,False,True,False,True,False]
    xlabels=["Max Air Temperature Bias","Max Air Temperature RMS error",
             "Min Air Temperature Bias","Min Air Temperature RMS error",
             "Precipitation Bias","Precipitation RMS error"]
    fignames=["tasmax_bias_histogram.png","tasmax_rms_histogram.png",
              "tasmin_bias_histogram.png","tasmin_rms_histogram.png",
              "pr_bias_histogram.png","pr_rms_histogram.png"]
    xlims=[[-2.5,2.5],[0,10],[-2.5,2.5],[0,10],[-400,400],[0,200]]
    for f,pz,xlab,xlim,fign in zip(filesearches,plotzeros,xlabels,xlims,fignames):
        plot_for(f,varname="data",plotzero=pz,trailing="n",
                     xlim=xlim,xlabel=xlab,title=None,figname=fign)
        
def prbias():
    starts=["CA_ncep_","SD_ncep_","SDmon_ncep_",
            "CA_narr_","SD_narr_","SDmon_narr_"]
    trail="pr_bias_map.nc"
    files=[f+trail for f in starts]
    
    plot_for(files=files,plotzero=True,xlim=[-350,350],varname="data",
             xlabel="Precip Bias (mm)",figname="pr_bias_histogram.png")
    

def wetfrac12():
    starts=["ncep_CA","ncep_SD","ncep_SDmon","obs_maurer",
            "narr_CA","narr_SD","narr_SDmon"]
    trail="_12km_pr_wet_frac.nc"
    files=[f+trail for f in starts]
    
    plot_for("*12km_pr_wet_frac.nc",files=files,plotzero=False,xlim=[0.2,1.0],varname="wetfrac",
             xlabel="Fraction Wet Days",figname="histo_wet_frac12.png",trailing="12",nbins=25,ylim=[0,20000])


def wetfrac():
    starts=["ncep_CA","ncep_SD","ncep_SDmon","obs_uw",
            "narr_CA","narr_SD","narr_SDmon","narr_SARe0"]
    trail="_4km_pr_wet_frac.nc"
    files=[f+trail for f in starts]
    
    plot_for("*4km_pr_wet_frac.nc",files=files,plotzero=False,xlim=[0.2,1.0],varname="wetfrac",
             xlabel="Fraction Wet Days",figname="histo_wet_frac.png",trailing="4",nbins=25)

def percentile(value):
    plot_for("*4km_pr_"+value+"percentile.nc",plotzero=False,xlim=[0,50],varname="precip",
             xlabel="99th Percentile Precip (mm)",figname="histo_99percentiles.png",
             trailing="4",nbins=50)

def ratiopercentile(value):
    
    search1="ncep_*4km_pr_"+value+"percentile.nc"
    search2="narr_*4km_pr_"+value+"percentile.nc"
    trailing="4"
    files=[f for f in glob.glob(search1) if (not re.match(".*dry_"+trailing+".*",f) 
                                        and not re.match(".*wet_"+trailing+".*",f) 
                                        and not re.match(".*hot_"+trailing+".*",f) 
                                        and not re.match(".*cold_"+trailing+".*",f))]
    files.extend([f for f in glob.glob(search2) if (not re.match(".*dry_"+trailing+".*",f) 
                                        and not re.match(".*wet_"+trailing+".*",f) 
                                        and not re.match(".*hot_"+trailing+".*",f) 
                                        and not re.match(".*cold_"+trailing+".*",f))])

    
    # files=glob.glob("CA*4km_pr_"+value+"percentile.nc")
    # files.extend(glob.glob("SD*4km_pr_"+value+"percentile.nc"))
    data=[]
    obs=mygis.read_nc("obs_uw_4km_pr_"+value+"percentile.nc","precip").data
    for f in files:
        data.append((mygis.read_nc(f,"precip").data-obs)/obs)
    
    plot_for(files=files,data=data,plotzero=True,xlim=[-1,1],varname="precip",
             xlabel="Fractional Error in "+value+"th Percentile Precip",
             figname="histo_"+value+"percentile_ratio.png",
             trailing="4",nbins=50)

        
def main(test):
    if test=="stat_v_obs":
        stat_v_obs()
    if re.match(".+percentile.*",test):
        percentile(test[:2])
    if re.match("percentile.*",test):
        percentile("99")
    if re.match(".+frac_percent.*",test):
        ratiopercentile(test[:2])
    if re.match("frac_percent.*",test):
        ratiopercentile("99")
    if test=="dryspell":
        dryspell()
    if test[:7]=="wetfrac":
        wetfrac12()
    if test=="prbias":
        prbias()
    

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='This is a template file for Python scripts. ')
        parser.add_argument('filename',action='store')
        parser.add_argument('-v', '--version',action='version',
                version='Template Parser 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()
    
        exit_code = main(args.filename)
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
