#!/usr/bin/env python
from matplotlib.pyplot import plot,ylabel,legend,clf,figure,xlim,ylim,yscale,xscale
from numpy import zeros
import numpy as np
import load_data, date_fun
import sys

useredgrey=True

prefix='chimney2d_realveg'
green=load_data.cols(prefix+'_basin.daily')
red=load_data.cols(prefix+'_upveg_basin.daily')
grey=load_data.cols('chimney2d_allpine_basin.daily')

greenmjd=date_fun.date2mjd(green[:,2],green[:,1],green[:,0],zeros(len(green[:,0])),zeros(len(green[:,0])),zeros(len(green[:,0])))
greendates=date_fun.mjd2datetime(greenmjd);
cur=red
mjd=date_fun.date2mjd(cur[:,2],cur[:,1],cur[:,0],zeros(len(cur[:,0])),zeros(len(cur[:,0])),zeros(len(cur[:,0])))
reddates=date_fun.mjd2datetime(mjd);
cur=grey
mjd=date_fun.date2mjd(cur[:,2],cur[:,1],cur[:,0],zeros(len(cur[:,0])),zeros(len(cur[:,0])),zeros(len(cur[:,0])))
greydates=date_fun.mjd2datetime(mjd);

curfig=figure()

label1='Real'
label2='Real*3'
label3='Homog'

##############################################
# PLOT Snow Water Equivalent
##############################################
plot(greendates,green[:,14],color='g',label=label1);
plot(reddates,red[:,14],color='r',label=label2);
plot(greydates,grey[:,14],color='grey',label=label3);
legend(loc=2);
ylabel('SWE (mm)');
curfig.autofmt_xdate()
xlim(reddates[0],reddates[-1]);
curfig.savefig(prefix+'_SWE_full.pdf')
xlim(reddates[-365],reddates[-1]);
curfig.savefig(prefix+'_SWE_subset.pdf')

##############################################
# PLOT Soil Moisture
##############################################
clf();
plot(greendates,0.5-green[:,7]/1000,color='g',label='Real');
plot(reddates,0.5-red[:,7]/1000,color='r',label='Real*3');
plot(greydates,0.5-grey[:,7]/1000,color='grey',label='Homog');
ylabel("Soil Moisture");
legend(loc=2);
curfig.autofmt_xdate()
xlim(reddates[0],reddates[-1]);
curfig.savefig(prefix+'_SMC_full.pdf')

xlim(reddates[-365],reddates[-1]);
curfig.savefig(prefix+'_SMC_subset.pdf')


##############################################
# PLOT Accumulated Transpiration
##############################################
clf();
accvar=0
plot(greendates,green[:,32],color='g',label=label1);
plot(reddates,red[:,32]+accvar,color='r',label=label2);
plot(greydates,grey[:,32]+accvar,color='grey',label=label3);
ylabel("Acc. Transp. (mm)");
legend(loc=2);
curfig.autofmt_xdate()
xlim(reddates[0],reddates[-1]);
yr=ylim()
ylim(green[-len(red[:,0]),32],yr[1])
curfig.savefig(prefix+'_transp_full.pdf')

xlim(reddates[-365],reddates[-1]);
curfig.savefig(prefix+'_transp_subset.pdf')


##############################################
# PLOT Sublimation
##############################################
clf();
accvar=0
greensub=green[:,30].cumsum()
redsub=red[:,30].cumsum()
greysub=grey[:,30].cumsum()
plot(greendates,greensub,color='g',label=label1);
plot(reddates,redsub+accvar,color='r',label=label2);
plot(greydates,greysub+accvar,color='grey',label=label3);
ylabel("Acc. Sublimation (mm)");
legend(loc=2);
curfig.autofmt_xdate()
xlim(reddates[0],reddates[-1]);
yr=ylim()
ylim(greensub[-len(redsub)],yr[1])
curfig.savefig(prefix+'_sublim_full.pdf')

xlim(reddates[-365],reddates[-1]);
curfig.savefig(prefix+'_sublim_subset.pdf')


##############################################
# PLOT Stream flow
##############################################
clf();
plot(greendates,green[:,16]+green[:,17],'-',color='g',label=label1);
plot(reddates,red[:,16]+red[:,17],color='r',label=label2);
plot(greydates,grey[:,16]+grey[:,17],color='grey',label=label3);
legend(loc=2);
ylabel("Streamflow (mm)");
curfig.autofmt_xdate()
xlim(reddates[0],reddates[-1]);
yscale('log')
curfig.savefig(prefix+'_Q_full.pdf')

xlim(reddates[-365],reddates[-1]);
curfig.savefig(prefix+'_Q_subset.pdf')

