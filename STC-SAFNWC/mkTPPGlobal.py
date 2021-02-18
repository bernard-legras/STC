#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the cold point temperature and pressure and store them in the TPP directory
as hdf5 files.

Created on Mon June 4 2018

@author: Bernard Legras
"""

import os
from ECMWF_N import ECMWF
from datetime import datetime, timedelta
import flammkuchen as fl
maindir = '/data/legras/flexpart_in/STC/ERA5/TPP'
#maindir = '/home/legras/sandbox/TPP'

date = datetime(2017,8,11,12)

while date < datetime(2017,8,11,18):
    print('processing ',date)
    outfile = date.strftime('TPP%y%m%d%H.hdf5')
    fullname = os.path.join(maindir,date.strftime('%Y/%m'),outfile)
    fdd = ECMWF('FULL-EA',date)
    fdd._get_T()
    fdd._mkp()
    fdd.close()
    fde = fdd.shift2west(-20)
    fdf = fde.extract(lonRange=[-10,160],latRange=[0,50],varss='All')
    fdf._CPT()
    fdf._WMO()
    tpp = {}
    tpp['Twmo'] = fdf.d2d['Twmo']
    tpp['pwmo'] = fdf.d2d['pwmo']
    tpp['Tcold'] = fdf.d2d['Tcold']
    tpp['pcold'] = fdf.d2d['pcold']
    tpp['nlon'] = fdf.nlon
    tpp['nlat'] = fdf.nlat
    tpp['lats'] = fdf.attr['lats']
    tpp['lons'] = fdf.attr['lons']
    tpp['date'] = fdd.date
    fl.save(fullname,tpp,compression='zlib')
    date += timedelta(hours=3)
