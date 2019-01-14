#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 29 May 2017

This script tests and demonstrates the patching features of geosat.py. 

@author: Bernard Legras
"""
import geosat
#import pickle,gzip
#import numpy as np
from datetime import datetime
#import matplotlib.pyplot as plt
#%%
date=datetime(year=2017,month=3,day=15,hour=19)
# Himawari read
hima=geosat.Himawari(date)
# read IR0
hima._get_IR0()
hima._get_var('WV_062')
hima._get_var('WV_073')
hima.close()
# MSG1 read
msg1=geosat.MSG1(date)
msg1._get_IR0()
msg1._get_var('WV_062')
msg1._get_var('WV_073')
msg1.close()
# grid projection and match
for grid in ['KTM','MesoInd','FullAMA']:
    gg = geosat.GeoGrid(grid)
    # associate grid and sat
    phima = geosat.SatGrid(hima,gg)
    pmsg1 = geosat.SatGrid(msg1,gg)
    # project by closest neigbour
    # and do angular correction
    for var in ['IR0','WV_062','WV_073']:
        phima._sat_togrid(var)
        pmsg1._sat_togrid(var)
    # angular correction for 'IR0'
    phima._sza_correc()
    pmsg1._sza_correc()
   
    # patching
    ppt = pmsg1.patch(phima,90.75,['IR0','WV_062','WV_073'])
    ppt.chart('IR0',txt='brightness temperature infrared 10.8')
    ppt.chart('WV_062',clim=[190,250],txt='6.2 µm')
    ppt.chart('WV_073',clim=[190,270],txt='7.3 µm')
