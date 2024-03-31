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
version = 'v2018.1-HVR'


# Get handle on satellite images for this date
# Standard GEO files
ah = geosat.MSG2(date)
# %%SAFNWC products
acth = SAFNWC_CT(date,'msg2')
achh = SAFNWC_CTTH(date,'msg2')
#acth2 = SAFNWC_CT(date,'himawari',BBname='ACCLIPBox',version=version)
#achh2 = SAFNWC_CTTH(date,'himawari',BBname='ACCLIPBox',version=version)

# %%Get handle on the ACCLIP grid
gg = geosat.GeoGrid('Kochi')
#gg2 = geosat.GeoGrid('ACCLIP_SAFBox')

# Associate satellites to grid
ph = geosat.SatGrid(ah,gg)
#ph2 = geosat.SatGrid(acth2,gg2)

# Read IR0 field
ah._get_IR0()
ah.close()

# Read CT field
acth._CT()
# Read PRESS and TEMP
achh._CTTH_PRESS()
achh._CTTH_TEMPER()
achh._CTTH_HEIGHT()
#
# Merge
ah._merge(acth)
ah._merge(achh)
del achh


# Project onto grid
varlist = ['IR0','CT','CTTH_PRESS','CTTH_TEMPER','CTTH_HEIGHT']
for var in  varlist:
   ph._sat_togrid(var)

#%%

# Plot images
ph.chart('IR0',txt='IR0 msg2')
ph.chart('CT',clim=(0,15),txt='CT msg2')
ph.chart('CTTH_PRESS',clim=(50,500),txt='CTTH_PRESS msg2')
ph.chart('CTTH_HEIGHT',clim=(0,20000),txt='CTTH_HEIGHT msg2')
ph.chart('CTTH_TEMPER',clim=(190,300),txt='CTTH_TEMPER msg2')


#%%
subgg = gg.subgrid([73,79,7,13])
ph.chart('IR0',txt='subgrid IR0',subgrid=subgg)
ph.chart('CT',clim=(0,15),txt='CT msg2',subgrid=subgg)
ph.chart('CTTH_HEIGHT',clim=(0,20000),txt='CTTH_HEIGHT msg2',subgrid=subgg)