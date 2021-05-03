#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the cold point temperature and pressure and store them in the TPP directory
as hdf5 files.

Created on Mon June 4 2018

@author: Bernard Legras
"""

import os
import numpy as np
from ECMWF_N import ECMWF
from datetime import datetime, timedelta
import argparse
import flammkuchen as fl
maindir = '/data/legras/flexpart_in/STC/ERA5/SAFNWP/TPP/LR'
#maindir = '/home/legras/sandbox/TPP'

parser = argparse.ArgumentParser()
parser.add_argument("-y","--year",type=int,help="year")
parser.add_argument("-m1","--month1",type=int,choices=1+np.arange(12),help="start month")
parser.add_argument("-d1","--day1",type=int,choices=1+np.arange(31),help="start day")
parser.add_argument("-m2","--month2",type=int,choices=1+np.arange(12),help="end month")
parser.add_argument("-d2","--day2",type=int,choices=1+np.arange(31),help="end day")
parser.add_argument("-h1","--hour1",type=int,choices=[0,3,6,9,12,15,18,21],help="start hour")
parser.add_argument("-h2","--hour2",type=int,choices=[0,3,6,9,12,15,18,21],help="end hour")
#parser.add_argument("-q","--quiet",type=str,choices=["y","n"],help= "quiet ys (y) or not (n)")

year = 2017
month1 = 8
day1 = 1
hour1 = 0
month2 = 9
day2 = 1
hour2 = 0

print('parsing arguments')
args = parser.parse_args()
if args.year is not None: year = args.year
if args.month1 is not None: month1 = args.month1
if args.month2 is not None: month2 = args.month2
if args.day1 is not None: day1 = args.day1
if args.day2 is not None: day2 = args.day2
if args.hour1 is not None: hour1 = args.hour1
if args.hour2 is not None: hour2 = args.hour2

date = datetime(year,month1,day1,hour1)

while date < datetime(year,month2,day2,hour2):
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
