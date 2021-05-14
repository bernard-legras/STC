#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on Thu Apr 23 00:28:34 2020

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime
#import pickle,gzip
#import matplotlib.pyplot as plt
from os.path import join

from ECMWF_N import ECMWF

date = datetime(2017,8,23,3)
dat = ECMWF('FULL-EA',date,exp='VOZ')

#date = datetime(2020,1,23,6)
#dat = ECMWF('OPZ',date)

#date = datetime(2017,8,10,0)
#dat = ECMWF('STC',date)

#%%
dat._get_T()
dat._get_var('VO')
dat._get_U()
dat._get_V()
dat._mkp()
dat._mkthet()
#%%
print('now call mkpv')
dat._mkpv()

#%%
#dats = dat.interpolPT([350,370,395,430],varList=['PV'],lonRange=(50,120),latRange=(0,50))
dats = dat.extract(varss='All',lonRange=(240,330),latRange=(-80,-30))
