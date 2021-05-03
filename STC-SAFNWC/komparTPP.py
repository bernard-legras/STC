#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script performs a series of tests on the WMO tropopause as
calculated from the ERA5

Created on Fri Feb 19 23:09:43 2021

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
from os.path import join
import flammkuchen as fl
from scipy.interpolate import RectBivariateSpline

def main():

    p0 = 100000
    H = 7.4
    mainHR = '/data/legras/flexpart_in/STC/ERA5/SAFNWP/TPP/HR'
    mainLR = '/data/legras/flexpart_in/STC/ERA5/SAFNWP/TPP/LR'

    """
    1 First test
    The high resolution WMO tropopause is comapared to the cold point tropopause and to the
    interpolated tropopause from the low reolution version at 1°x1°

    It is found that the high resoluton version exhibits a small number of spurious low level
    tropopauses in the monsoon region which are not obtainde in the llow resolution version.

    """
    date0 = datetime(2017,8,10,12)
    date0 = datetime(2017,8,31,0)
    date0 = datetime(2017,8,6,0)
    date0 = datetime(2017,8,13,3)

    date = date0
    file = date.strftime('TPP%y%m%d%H.hdf5')
    # get high-resolution version
    HRname = join(mainHR,date.strftime('%Y/%m'),file)
    tppHR = fl.load(HRname)
    print(np.min(tppHR['Tcold']),np.max(tppHR['Tcold']))
    print(np.min(tppHR['pcold']),np.max(tppHR['pcold']))
    tppHR['zcold'] = - H * np.log(tppHR['pcold']/p0)
    tppHR['zwmo'] = - H * np.log(tppHR['pwmo']/p0)
    # get low resolution version
    LRname = join(mainLR,date.strftime('%Y/%m'),file)
    tppLR = fl.load(LRname)
    tppLR['zcold'] = - H * np.log(tppLR['pcold']/p0)
    tppLR['zwmo'] = - H * np.log(tppLR['pwmo']/p0)
    # Interpolate from lR to HR
    tppZI = RectBivariateSpline(tppLR['lats'],tppLR['lons'],tppLR['zwmo'],kx=3,ky=3)
    tppTI = RectBivariateSpline(tppLR['lats'],tppLR['lons'],tppLR['Twmo'],kx=3,ky=3)
    tppHR['zwmo_smooth'] = tppZI(tppHR['lats'],tppHR['lons'])
    tppHR['Twmo_smooth'] = tppTI(tppHR['lats'],tppHR['lons'])

    # Differences HR - LR and WMO to cold point
    tppHR['dz'] = tppHR['zwmo'] - tppHR['zwmo_smooth']
    tppHR['dzc'] = tppHR['zwmo'] - tppHR['zcold']

    show(tppHR,'zcold',txt='Z cold point tropopause (km) HR')
    show(tppHR,'zwmo',txt='Z WMO tropopause (km) HR')

    show(tppLR,'zcold',txt='Z cold point tropopause (km) LR')
    show(tppLR,'zwmo',txt='Z WMO tropopause (km) LR')
    show(tppHR,'zwmo_smooth',txt='Z WMO tropopause (km) smoother')
    show(tppHR,'dzc',txt='Delta Z WMO - cold HR (km)')
    show(tppHR,'dz',txt='Delta Z WMO HR - smoothed LR (km)')
    #plt.hist(tppHR['zcold']);plt.show()
    #plt.hist(tppHR['zwmo']);plt.show()
    #plt.hist(tppHR['zwmo_smooth']);plt.show()

    """
    2 Second test
    Statistics on the anomalies in the 0-32° region
    """
    date0 = datetime(2017,8,1,0)
    date = date0
    nLR = []
    nHR = []
    while date < datetime(2017,9,1,0):
        file = date.strftime('TPP%y%m%d%H.hdf5')
        HRname = join(mainHR,date.strftime('%Y/%m'),file)
        tppHR = fl.load(HRname)
        tppHR['zwmo'] = - H * np.log(tppHR['pwmo']/p0)
        LRname = join(mainLR,date.strftime('%Y/%m'),file)
        tppLR = fl.load(LRname)
        tppLR['zwmo'] = - H * np.log(tppLR['pwmo']/p0)
        tppZI = RectBivariateSpline(tppLR['lats'],tppLR['lons'],tppLR['zwmo'],kx=3,ky=3)
        tppHR['zwmo_smooth'] = tppZI(tppHR['lats'],tppHR['lons'])
        nHR.append(np.sum(tppHR['zwmo'][:128,:] < 12))
        nLR.append(np.sum(tppHR['zwmo_smooth'][:128,:] < 12))
        print(date.strftime('%d %B %H UTC'),nLR[-1],nHR[-1],np.min(tppLR['pwmo'])/100)
        date += timedelta(hours=3)
    return([nLR,nHR])

def show(dic,var,txt=None):
    fig = plt.figure(figsize=(12,3))
    ax = plt.axes()
    im = ax.imshow(dic[var],extent=[-10,160,0,50],origin='lower',cmap ='jet')
    if txt is not None: ax.set_title(txt)
    plt.colorbar(im)
    plt.show()

if __name__ == '__main__':
    main()
