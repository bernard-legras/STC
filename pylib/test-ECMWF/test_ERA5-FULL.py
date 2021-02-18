#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 22:45:57 2018

@author: Bernard Legras
"""
from ECMWF_N import ECMWF
from datetime import datetime

date = datetime(2017,8,11,12)

dat = ECMWF('FULL-EA',date)
dat._get_T()
dat._mkp()
dat._mkz()
#dat._CPT()
#dat.show('zcold')
[lapse,dz]=dat._WMO()
dat.show('zwmo')
