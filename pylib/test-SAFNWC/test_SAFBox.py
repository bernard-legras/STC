#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is testing the SAFBox

Created on Fri Jun  8 13:33:15 2018

@author: Bernard Legras
"""

import geosat
from SAFNWCnc import SAFNWC_CTTH
from datetime import datetime

#date = datetime(2017,7,22,3,45)
date = datetime(2017,8,8,9,00)
dat = SAFNWC_CTTH(date,'msg1',BBname='SAFBox',version='v2018.1-HVR')
gg = geosat.GeoGrid('FullAMA_SAFBox')
dat._CTTH_PRESS()
p1 = geosat.SatGrid(dat,gg)
p1._sat_togrid('CTTH_PRESS')
p1.chart('CTTH_PRESS')
