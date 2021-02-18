#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the cold point temperature and pressure and store them in the TPP directory
as hdf5 files.

Created on Mon June 4 2018
Modified on Mon 15 Februray 2020 to add WMO topopause

@author: Bernard Legras
"""

import os
from ECMWF_N import ECMWF
from datetime import datetime, timedelta
import flammkuchen as fl

maindir = '/data/legras/flexpart_in/STC/ERA5/TPP/HR'

date = datetime(2017,8,1,0)

while date < datetime(2017,9,1,0):
    print('processing ',date)
    outfile = date.strftime('TPP%y%m%d%H.hdf5')
    fullname = os.path.join(maindir,date.strftime('%Y/%m'),outfile)
    fdd = ECMWF('STC',date)
    fdd._get_T()
    fdd._mkp()
    fdd._CPT()
    fdd._WMO()
    fdd.close()
    tpp = {}
    tpp['Twmo'] = fdd.d2d['Twmo']
    tpp['pwmo'] = fdd.d2d['pwmo']
    tpp['Tcold'] = fdd.d2d['Tcold']
    tpp['pcold'] = fdd.d2d['pcold']
    tpp['nlon'] = fdd.nlon
    tpp['nlat'] = fdd.nlat
    tpp['lats'] = fdd.attr['lats']
    tpp['lons'] = fdd.attr['lons']
    tpp['date'] = fdd.date
    fl.save(fullname,tpp,compression='zlib')
    date += timedelta(hours=1)
