#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the cold point temperature and pressure and store them in the TPP directory
as hdf5 files.

Created on Mon June 4 2018

@author: Bernard Legras
"""

import os
import shutil
import pygrib
from datetime import datetime, timedelta
import deepdish as dd

ERA5 ='/dsk2/ERA5/STC'
TPPdir = os.path.join(ERA5,'TPP')
ENPsourceDir = os.path.join(ERA5,'ENPsource')
ENPdir = os.path.join(ERA5,'ENP')
tmpFile = 'tpp.grib'

date = datetime(2017,9,18)
# back one increment to make it at the beginning of the loop
date -= timedelta(hours=1)

while date < datetime(2017,10,1):
    #increment the date
    date += timedelta(hours=1)
    print('processing ',date)

    # read the tropopause data
    tppFile = date.strftime('TPP%y%m%d%H.hdf5')
    TPPfullName = os.path.join(TPPdir,date.strftime('%Y/%m'),tppFile)
    try:
        tpp = dd.io.load(TPPfullName)
    except:
        print('tpp file not read')
        print(TPPfullName)
        continue

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
        mmm['values'] = tpp['Tcold']
        mmm.__setitem__('paramId',21)
        msg = mmm.tostring()
        grbOut.write(msg)

        # build and write the second message
        mmm['values'] = tpp['pcold']/100
        mmm.__setitem__('paramId',22)
        msg = mmm.tostring()
        grbOut.write(msg)

    # copy the ENP file and add tropopause messages
    shutil.copy(ENPsourceFile,ENPtargetFile)
    with open(ENPtargetFile,'ab') as target, open(tmpFile,'rb') as tppFile:
        target.write(tppFile.read())
