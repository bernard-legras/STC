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
#import matplotlib.pyplot as plt
#%%
date=datetime(year=2022,month=1,day=27,hour=15)
# Himawari read and plot
ah=geosat.Himawari(date)
gg=geosat.GeoGrid('FullAMA')
ah._get_IR0()
ph=geosat.SatGrid(ah,gg)
ph._sat_togrid('IR0')
ah.show('IR0')
ph.chart('IR0',txt='BT himawari')
# MSG1 read and plot
a1=geosat.MSG1(date)
a1._get_IR0()
p1=geosat.SatGrid(a1,gg)
p1._sat_togrid('IR0')
a1.show('IR0')
p1.chart('IR0',txt='BT msg1')
# MSG3 read and plot
a3=geosat.MSG0(date)
a3._get_IR0()
p3=geosat.SatGrid(a3,gg)
p3._sat_togrid('IR0')
a3.show('IR0')
p3.chart('IR0',txt='BT msg3')
# difference and plot
dd=geosat.GridField(gg)
dd.var['IR0']=p1.var['IR0']-ph.var['IR0']
dd.chart('IR0',clim=(-10,10),txt='BT msg1-himawari')
dd.var['IR0']=p3.var['IR0']-p1.var['IR0']
dd.chart('IR0',clim=(-10,10),txt='BT msg3-msg1')

##%% Comparison of the image sizes
## load coordinate mask
#lonlat3=pickle.load(gzip.open('msg3/lonlat.pkl','rb'))
## load MSG3 image
#a3=geosat.MSG3(date)
#a3._get_IR0()
## compares size of images as a function of line
#diff13=np.empty(3712,dtype=int)
#diff33=np.empty(3712,dtype=int)
#for i in range(3712):
#    diff13[i] = len(a1.var['IR0'][i,:].compressed())-len(a3.var['IR0'][i,:].compressed())
#    diff33[i] = len(lonlat3['lon'][i,:].compressed())-len(a3.var['IR0'][i,:].compressed())
#plt.plot(diff13)
#plt.show()
#plt.plot(diff33)
#plt.show()
#%%
# application of the zenit angle correction and plot
ph._sza_correc()
p1._sza_correc()
p3._sza_correc()
ph.chart('IR0',txt='BT himawari corrected')
p1.chart('IR0',txt='BT msg1 corrected')
p3.chart('IR0',txt='BT msg0 corrected')
# difference and plot
dd.var['IR0']=p1.var['IR0']-ph.var['IR0']
dd.chart('IR0',clim=(-10,10),txt='BT msg1-himawari corrected')
dd.var['IR0']=p3.var['IR0']-p1.var['IR0']
dd.chart('IR0',clim=(-10,10),txt='BT msg0-msg1 corrected')
#%% Test patching
patched=p1.patch(ph,90,'IR0')
patched.chart('IR0',txt='patched image')
