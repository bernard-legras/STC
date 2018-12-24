#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script tests the behavior of geosat with a MSG1 image with missing lines

Created on Sat Jun 10 16:16:30 2017

@author: Bernard Legras
"""
import geosat
from datetime import datetime

date = datetime(2017,6,8,7)

ah = geosat.MSG1(date)
gg = geosat.GeoGrid('MesoInd')

ah._get_IR0()
ph=geosat.SatGrid(ah,gg)
ph._sat_togrid('IR0')
ah.show('IR0')
ph.chart('IR0',txt='BT MSG1')