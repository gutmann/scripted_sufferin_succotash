#!/usr/bin/env python
# encoding: utf-8
"""
smc_tsoil_corr.py

Created by Ethan Gutmann on 2011-05-06.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""

import numpy as np
from ols import ols
import matplotlib.pyplot as plt

def calc(smc,ts,plot=None,force=None):
	n=smc.size
	dummy=np.arange(n)
	
	smc_init=ols(smc,dummy)
	smctest=smc-(smc_init.b[0]+smc_init.b[1]*dummy)
	
	ts_init=ols(ts,dummy)
	tstest=ts-(ts_init.b[0]+ts_init.b[1]*dummy)
	
	smc_ts=ols(smctest,tstest)
	
	if plot!=None:
		plt.clf()
		plt.plot(tstest/1000)
		plt.plot(smctest)
		if force==None:
			plt.plot(smctest-(tstest*smc_ts.b[1]))
		else:
			plt.plot(smctest-(tstest*force))
		print(smc_ts.R2)
	return smc_ts.b
