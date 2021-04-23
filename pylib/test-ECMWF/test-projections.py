#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on Sun May 31 17:35:02 2020

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from ECMWF_N import ECMWF

dat = ECMWF('FULL-EA',datetime(2020,1,21,6),exp='VOZ')
dat._get_var('VO')
dat._get_var('T')

#%% Example of a plot that is contained within the southern west hemisphere

dats = dat.extract(varss=['T','VO'],lonRange=(240,300),latRange=(-80,-30))

dats.show('VO',49)
dats.show('VO',49,projec='ortho',xylim=True)

# Same plot with a projection where the shift is -360 everywhere (proj extent, boundaries)
ortho =ccrs.Orthographic(270-360,-55)
extent = (240-360,300-360,-80,-30)
ax = plt.axes(projection=ortho)
ax.imshow(dats.var['VO'][49,...],transform=ccrs.PlateCarree(),
                  extent=extent, origin='lower',interpolation='nearest',
                  cmap='jet')
x1,_ = ortho.transform_point(240-360,-30,ccrs.Geodetic())
_,y1 = ortho.transform_point(240-360,-80,ccrs.Geodetic())
x2,_ = ortho.transform_point(300-360,-30,ccrs.Geodetic())
_,y2 = ortho.transform_point(270-360,-30,ccrs.Geodetic())
ax.gridlines()
ax.coastlines()
ax.set_xlim(x1,x2)
ax.set_ylim(y1,y2)
plt.show()

#%% Example of a plot that is contained within the southern east hemisphere

dats = dat.extract(varss=['T','VO'],lonRange=(30,150),latRange=(-80,-30))

dats.show('VO',49)
dats.show('VO',49,projec='ortho',xylim=True)

# Same plot with a projection where the shift is -360 everywhere (proj extent, boundaries)
ortho =ccrs.Orthographic(90,-55)
extent = (30,150,-80,-30)
ax = plt.axes(projection=ortho)
ax.imshow(dats.var['VO'][49,...],transform=ccrs.PlateCarree(),
                  extent=extent, origin='lower',interpolation='nearest',
                  cmap='jet')
x1,_ = ortho.transform_point(30,-30,ccrs.Geodetic())
_,y1 = ortho.transform_point(30,-80,ccrs.Geodetic())
x2,_ = ortho.transform_point(150,-30,ccrs.Geodetic())
_,y2 = ortho.transform_point(90,-30,ccrs.Geodetic())
ax.gridlines()
ax.coastlines()
ax.set_xlim(x1,x2)
ax.set_ylim(y1,y2)
plt.show()

#%% Example of a plot that straddles the date line in the southern hemisphere

dats = dat.extract(varss=['T','VO'],lonRange=(120,240),latRange=(-80,-30))

dats.show('VO',49)
dats.show('VO',49,projec='ortho',xylim=True)

#%% Example of a plot that straddles the date line in the northern hemisphere

dats = dat.extract(varss=['T','VO'],lonRange=(120,240),latRange=(30,80))

dats.show('VO',49)
dats.show('VO',49,projec='ortho',xylim=True)
# other projections
dats.show('VO',49,projec='lambert',txt='Lambert Conformal',xylim=True)
dats.show('VO',49,projec='azimuthalequi',txt='Azimuthal Equidistant',xylim=True)
dats.show('VO',49,projec='azimuthalequi',txt='Azimuthal Equidistant',xylim=True)
dats.show('VO',49,projec='nearside',txt='Nearside Perspective',xylim=True)

#%% Example of a plot over 180° in longitude

dats = dat.extract(varss=['T','VO'],lonRange=(60,240),latRange=(30,80))

dats.show('VO',49)
dats.show('VO',49,projec='ortho',xylim=True)

#%% Example of a plot over more than 180° in longitude

dats = dat.extract(varss=['T','VO'],lonRange=(30,270),latRange=(30,80))

dats.show('VO',49)
# In this case xylim produces truncature or even errors if some points are hidden on the other side
dats.show('VO',49,projec='ortho',xylim=True)
# It is then better to set the projection above the pole
dats.show('VO',49,projec='ortho',polar=True)
# Still better to use xylim for lambert
dats.show('VO',49,projec='lambert',xylim=True)