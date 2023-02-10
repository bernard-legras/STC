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
#%%geosat.Himawari(date
date=datetime(year=2022,month=1,day=27,hour=15)

age = geosat.GOESE(date)
ge = geosat.GeoGrid('GOESEFull')
age._get_IR0()
pge = geosat.SatGrid(age,ge)
pge._sat_togrid('IR0')
pge.chart('IR0')

agw = geosat.GOESW(date)
gw = geosat.GeoGrid('GOESWFull')
agw._get_IR0()
pgw = geosat.SatGrid(agw,gw)
pgw._sat_togrid('IR0')
pgw.chart('IR0')