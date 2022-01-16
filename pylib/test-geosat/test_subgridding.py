#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 19:03:20 2017

Demonstration and test of subgridding

@author: legras
"""
import geosat
#import pickle,gzip
#import numpy as np
from datetime import datetime
#import matplotlib.pyplot as plt
#%%
date=datetime(year=2017,month=3,day=15,hour=19)
# Himawari read and plot
ah=geosat.Himawari(date)
gg=geosat.GeoGrid('MesoInd')
ah._get_IR0()
ph=geosat.SatGrid(ah,gg)
ph._sat_togrid('IR0')
ah.show('IR0')
ph.chart('IR0',txt='BT himawari')

# exemple of subgridding and plot with this subgrid
subgg=gg.subgrid([73,97,11,35])
print(subgg.box_range)
print(subgg.shapexy)
ph.chart('IR0',txt='subgrid test',subgrid=subgg)