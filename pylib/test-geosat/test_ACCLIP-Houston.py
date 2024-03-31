#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 16:53:21 2023

@author: Bernard Legras
"""
import geosat
from SAFNWCnc import SAFNWC_CT, SAFNWC_CTTH
from datetime import datetime

# Date 
date = datetime(2022,8,31,15)

# Get handle on satellite images for this date
# Standard GEO files
ah = geosat.Himawari(date)
ag = geosat.GOESW(date)
am = geosat.GOESE(date)
# SAFNWC products
acth = SAFNWC_CT(date,'himawari')
actg = SAFNWC_CT(date,'goesw')
actm = SAFNWC_CT(date,'goese')
achh = SAFNWC_CTTH(date,'himawari')
achg = SAFNWC_CTTH(date,'goesw')
achm = SAFNWC_CTTH(date,'goese')

# Get handle on the ACCLIP grid
gg = geosat.GeoGrid('ACCLIP-Houston')

# Associate satellites to grid
ph = geosat.SatGrid(ah,gg)
pg = geosat.SatGrid(ag,gg)
pm = geosat.SatGrid(am,gg)

# Read IR0 field
ah._get_IR0()
ag._get_IR0()
am._get_IR0()
ah.close()
ag.close()
am.close()

# Read CT field
acth._CT()
actg._CT()
actm._CT()
# Merge
ah._merge(acth)
ag._merge(actg)
am._merge(actm)
del acth, actg, actm

# Read PRESS and TEMP
achh._CTTH_PRESS()
achg._CTTH_PRESS()
achm._CTTH_PRESS()
achh._CTTH_TEMPER()
achg._CTTH_TEMPER()
achm._CTTH_TEMPER()
# Merge
ah._merge(achh)
ag._merge(achg)
am._merge(achm)
del achh, achg, achm

# Project onto grid
varlist = ['IR0','CT','CTTH_PRESS','CTTH_TEMPER']
for var in  varlist:
    ph._sat_togrid(var)
    pg._sat_togrid(var)
    pm._sat_togrid(var)

# Generates composite by patching from west to east 
patch1 = ph.patch(pg,182.8,varlist)
patch = patch1.patch(pm,254,varlist)

# Plot composite
patch.chart('IR0',txt='IR0 composite')
patch.chart('CT',clim=(0,15),txt='CT composite')
patch.chart('CTTH_PRESS',clim=(50,500),txt='CTTH_PRESS composite')
patch.chart('CTTH_TEMPER',clim=(190,300),txt='CTTH_TEMPER composite')




