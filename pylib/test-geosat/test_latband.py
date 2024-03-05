#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:52:48 2017

This script tests and demonstrates the main features of geosat.py.

@author: Bernard Legras
"""
import geosat
#import pickle,gzip
#import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
#%%geosat.Himawari(date
date=datetime(year=2022,month=1,day=27,hour=15)

age = geosat.GOESE(date)
agw = geosat.GOESW(date)
am0 = geosat.MSG0(date)
am1 = geosat.MSG1(date)
ah = geosat.Himawari(date)
band = geosat.GeoGrid('LatBand1')
pge = geosat.SatGrid(age,band)
pgw = geosat.SatGrid(agw,band)
pm0 = geosat.SatGrid(am0,band)
pm1 = geosat.SatGrid(am1,band)
ph = geosat.SatGrid(ah,band)
age._get_IR0()
agw._get_IR0()
am0._get_IR0()
am1._get_IR0()
ah._get_IR0()
pge._sat_togrid('IR0')
pgw._sat_togrid('IR0')
pm0._sat_togrid('IR0')
pm1._sat_togrid('IR0')
ph._sat_togrid('IR0')
#%%
patch1 = pgw.patch(pge,-106,'IR0')
patch2 = patch1.patch(pm0,-38,'IR0')
patch3 = patch2.patch(pm1,21,'IR0')
patch4 = patch3.patch(ph,90,'IR0')
patchf = patch4.patch(pgw,181,'IR0')
#patchf = patch4
#%%
patchf.chart('IR0',txt='Full band IR0 composit')
plt.imshow(patchf.var['IR0'],extent=(-160,199.99,-35,10),origin='lower',cmap='jet',clim=(200,300))
