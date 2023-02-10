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
from datetime import datetime, timedelta
#import matplotlib.pyplot as plt
#%%geosat.Himawari(date
date=datetime(year=2022,month=1,day=27,hour=15)

try: age = geosat.GOESE(date)
except:
    try: age = geosat.GOESE(date+timedelta(hours=1))
    except: age = geosat.GOESE(date-timedelta(hours=1))
try: agw = geosat.GOESW(date)
except:
    try: agw = geosat.GOESW(date+timedelta(hours=1))
    except: agw = geosat.GOESW(date-timedelta(hours=1))
try: am0 = geosat.MSG9(date)
except:
    try: am0 = geosat.MSG0(date+timedelta(hours=1))
    except: am0 = geosat.MSG0(date-timedelta(hours=1))
try: am1 = geosat.MSG1(date)
except:
    try: am1 = geosat.MSG1(date+timedelta(hours=1))
    except: am1 = geosat.MSG1(date-timedelta(hours=1))
try: ah = geosat.Himawari(date)
except:
    try: ah = geosat.Himawari(date+timedelta(hours=1))
    except: ah = geosat.Himawari(date-timedelta(hours=1))

band = geosat.GeoGrid('LatBand1')
pge = geosat.SatGrid(age,band)
pgw = geosat.SatGrid(agw,band)
pm0 = geosat.SatGrid(am0,band)
pm1 = geosat.SatGrid(am1,band)
ph = geosat.SatGrid(ah,band)
agw._mk_Ash()
age._mk_Ash()
am0._mk_Ash()
am1._mk_Ash()
ah._mk_Ash()
pge._sat_togrid('Ash')
pgw._sat_togrid('Ash')
pm0._sat_togrid('Ash')
pm1._sat_togrid('Ash')
ph._sat_togrid('Ash')
#%%
patch1 = pgw.patch(pge,-106,'Ash')
patch2 = patch1.patch(pm0,-38,'Ash')
patch3 = patch2.patch(pm1,21,'Ash')
patchf = patch3.patch(ph,90,'Ash')
patchf.chart('Ash',txt='Full band Ash composit')