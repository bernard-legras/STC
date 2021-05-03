#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add the CPT and WMO tropopopause to the already built ENP file
containing the cold point tropopause but remove it first as it was badly oriented

The new temperature and pressure messages are added at the end of the file
with four parameter values
76, 77 for CPT tropopause
81, 82 for WMO tropopause

The latitudes are inverted to reverse the inversion made in the reading code
of ECMWF_N. It is expected that latituds are ranked from north to south in grib
files originating from ECMWF.

Created on Sat 20 Feb 2021

@author: Bernard Legras
"""

import os
import shutil
import numpy as np
import pygrib
from datetime import datetime, timedelta
from scipy.interpolate import RectBivariateSpline
import flammkuchen as fl
from subprocess import call

ERA5 ='/work_users/b.legras/ERA5/STC'
TPPdir = os.path.join(ERA5,'TPP/LR')
ENPsourceDir = os.path.join(ERA5,'ENPsource')
ENPdir = os.path.join(ERA5,'ENP')
tmpFile = 'tpp.grib'

date = datetime(2017,9,1,0)

if (date.hour%3) != 0:
    print('hour must be a multiple of 3')
    raise ValueError
else:
    # read first tropopause file
    tppFile = date.strftime('TPP%y%m%d%H.hdf5')
    TPPfullName = os.path.join(TPPdir,date.strftime('%Y/%m'),tppFile)
    tppN = fl.load(TPPfullName)

w1 = 1/3
w2 = 2/3
lons = np.arange(-10,160.1,0.25)
lats = np.arange(0,50.1,0.25)

while date < datetime(2017,10,1,0):

    print('processing ',date)

    # Read of the tropopause data

    if (date.hour%3) == 0:
        tppO = tppN.copy()
        # read the tropopause data for the next 3h slot
        tppFile = date.strftime('TPP%y%m%d%H.hdf5')
        TPPfullName = os.path.join(TPPdir,date.strftime('%Y/%m'),tppFile)
        tppN = fl.load(TPPfullName)
    # Time nterpolation
    if (date.hour%3) == 0:
        Twmo = tppO['Twmo']
        pwmo = tppO['pwmo']
        Tcpt = tppO['Tcold']
        pcpt = tppO['Tcold']
    elif (date.hour%3) == 1:
        Twmo = w2*tppO['Twmo'] + w1*tppN['Twmo']
        pwmo = w2*tppO['pwmo'] + w1*tppN['pwmo']
        Tcpt = w2*tppO['Tcold'] + w1*tppN['Tcold']
        pcpt = w2*tppO['pcold'] + w1*tppN['pcold']
    elif (date.hour%3) == 2:
        Twmo = w1*tppO['Twmo'] + w2*tppN['Twmo']
        pwmo = w1*tppO['pwmo'] + w2*tppN['pwmo']
        Tcpt = w1*tppO['Tcold'] + w2*tppN['Tcold']
        pcpt = w1*tppO['pcold'] + w2*tppN['pcold']

    # space interpolation
    tppTI = RectBivariateSpline(tppN['lats'],tppN['lons'],Twmo,kx=3,ky=3)
    TwmoHR = tppTI(lats,lons)
    tppZI  = RectBivariateSpline(tppN['lats'],tppN['lons'],np.log(pwmo),kx=3,ky=3)
    pwmoHR = np.exp(tppZI(lats,lons))/100
    tppTI = RectBivariateSpline(tppN['lats'],tppN['lons'],Tcpt,kx=3,ky=3)
    TcptHR = tppTI(lats,lons)
    tppZI  = RectBivariateSpline(tppN['lats'],tppN['lons'],np.log(pcpt),kx=3,ky=3)
    pcptHR = np.exp(tppZI(lats,lons))/100

    # Invert the latitudes for all four fields
    TwmoHR = TwmoHR[::-1,:]
    pwmoHR = pwmoHR[::-1,:]
    TcptHR = TcptHR[::-1,:]
    pcptHR = pcptHR[::-1,:]

    # define the source and target ENP files
    ENPfile = date.strftime('ENP%y%m%d%H')
    ENPsourceFile = os.path.join(ENPsourceDir,date.strftime('%Y/%m'),ENPfile)
    ENPtargetFile = os.path.join(ENPdir,date.strftime('%Y/%m'),ENPfile)

    #open the source file and extract the first message
    try:
        fid = pygrib.open(ENPsourceFile)
    except:
        print('ENP source file not found')
        continue
    mmm = fid.read(1)[0]
    fid.close()

    # open temp file for tropopause messages
    with open(tmpFile,'wb') as grbOut:

        # build and write first message
        mmm['values'] = TcptHR
        mmm.__setitem__('paramId',84)
        msg = mmm.tostring()
        grbOut.write(msg)

        # build and write the second message
        mmm['values'] = pcptHR
        mmm.__setitem__('paramId',85)
        msg = mmm.tostring()
        grbOut.write(msg)

        # build and write first message
        mmm['values'] = TwmoHR
        mmm.__setitem__('paramId',81)
        msg = mmm.tostring()
        grbOut.write(msg)

        # build and write the second message
        mmm['values'] = pwmoHR
        mmm.__setitem__('paramId',82)
        msg = mmm.tostring()
        grbOut.write(msg)

    # copy the surface data of the ENP file
    call(['grib_copy','-w','count=1/2/3/4/5/6/7/8',ENPsourceFile,'surf.grb'])
    # copy the pressure levels of the ENP file
    call(['grib_copy','-w','levtype=pl',ENPsourceFile,'pl.grb'])
    # make final copy
    shutil.copy('surf.grb',ENPtargetFile)
    with open(ENPtargetFile,'ab') as target, open(tmpFile,'rb') as tppFile:
        target.write(tppFile.read())
    with open(ENPtargetFile,'ab') as target, open('pl.grb','rb') as plFile:
        target.write(plFile.read())

    # increment the date"
    date += timedelta(hours=1)
