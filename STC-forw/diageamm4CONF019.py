#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plot of the mean properties of the heating rates and of the vertical velocities
during July august 2017, using the monthly archive of ERA5

Created on Sat Feb  9 04:02:06 2019

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
#import pickle,gzip
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from os.path import join
from eaMM import eaMM
import constants as cst

dmesh = cst.REarth*np.deg2rad(0.5)

#%%
# Color list with 20 colors      
listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']            
mymap=colors.ListedColormap(listcolors)

#%%
dat = eaMM(datetime(2017,7,1))
dat._getv('T')
dat._mkp()
dat._mkthet()
dat._mkz()
dat._getv('W')
dat._getv('D')
dat._getv('U')
dat._getv('V')
dat._getv('ASSWR')
dat._getv('ASLWR')
dat.var['AS'] = dat.var['ASSWR'] + dat.var['ASLWR']
dat8 = eaMM(datetime(2017,8,1))
dat8._getv('T')
dat8._mkp()
dat8._mkthet()
dat8._mkz()
dat8._getv('W')
dat8._getv('U')
dat8._getv('V')
dat8._getv('D')
dat8._getv('ASSWR')
dat8._getv('ASLWR')
dat8.var['AS'] = dat8.var['ASSWR'] + dat8.var['ASLWR']
# averaging July and August
for var in ['T','P','PT','W','D','AS','U','V','SP','Z']:
    dat.var[var] += dat8.var[var]
    dat.var[var] *= 0.5
del dat8
dat.var['ASPT'] = dat.var['AS']*(dat.var['P']/cst.p0)**(-cst.kappa)
dat.var['M'] = cst.Cp * dat.var['T'] + cst.g * dat.var['Z']

#%%
# colormesh wants edges as axes, do it in a quick and dirty way as we are
# only plotting a subset
# generate half coordinate mesh 
xedges = 0.5*(dat.lons[1:]+dat.lons[:-1])
jy = 240
print(dat.lats[jy])
ixr = (340,700)
ktr = (50,130)
tedges = 0.25*(dat.var['PT'][1:,jy,1:]+dat.var['PT'][:-1,jy,1:]
                 +dat.var['PT'][1:,jy,:-1]+dat.var['PT'][:-1,jy,:-1])
fs = 16
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(221)
im = ax.pcolormesh(xedges[ixr[0]-1:ixr[1]],tedges[ktr[0]-1:ktr[1],ixr[0]-1:ixr[1]],
                dat.var['ASPT'][ktr[0]:ktr[1],jy,ixr[0]:ixr[1]],cmap=mymap,clim=(-3,2))
ax.contour(np.broadcast_to(dat.lons[ixr[0]:ixr[1]],(ktr[1]-ktr[0],ixr[1]-ixr[0])),
            dat.var['PT'][ktr[0]:ktr[1],jy,ixr[0]:ixr[1]],
            dat.var['ASPT'][ktr[0]:ktr[1],jy,ixr[0]:ixr[1]],[0,])
plt.ylim(330,400)
plt.xlim(30,120)
plt.tick_params(labelsize=fs)
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=fs)
cbar.ax.set_title('K/day',size=fs)
plt.xlabel('longitude',fontsize=fs)
plt.ylabel('potential temperature (K)',fontsize=fs)
plt.title(r'Radiative D$\theta$/Dt at 30N',fontsize=fs)
#plt.show()

jy = 210
tedges = 0.25*(dat.var['PT'][1:,jy,1:]+dat.var['PT'][:-1,jy,1:]
                 +dat.var['PT'][1:,jy,:-1]+dat.var['PT'][:-1,jy,:-1])
ax = fig.add_subplot(222)
im=ax.contourf(np.broadcast_to(dat.lons[ixr[0]:ixr[1]],(ktr[1]-ktr[0],ixr[1]-ixr[0])),
             dat.var['PT'][ktr[0]:ktr[1],jy,ixr[0]:ixr[1]],
             dat.var['ASPT'][ktr[0]:ktr[1],jy,ixr[0]:ixr[1]],np.linspace(-3,2,11),
             cmap=mymap,extend='both')
ax.contour(np.broadcast_to(dat.lons[ixr[0]:ixr[1]],(ktr[1]-ktr[0],ixr[1]-ixr[0])),
             dat.var['PT'][ktr[0]:ktr[1],jy,ixr[0]:ixr[1]],
             dat.var['ASPT'][ktr[0]:ktr[1],jy,ixr[0]:ixr[1]],[0,])
plt.ylim(330,400)
plt.xlim(30,120)
plt.tick_params(labelsize=fs)
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=fs)
cbar.ax.set_title('K/day',size=fs)
plt.xlabel('longitude',fontsize=fs)
plt.ylabel('potential temperature (K)',fontsize=fs)
plt.title(r'Radiative D$\theta$/Dt at 15N',fontsize=fs)

# Vertical velocity as heating 

jy = 240
pspad = np.pad(dat.var['SP'][jy,:],(1,1),'wrap')
bl = 0.5*(dat.bl[1:]+dat.bl[:-1])
U = 0.5*(dat.var['U'][1:,jy,:]+dat.var['U'][:-1,jy,:])
V = 0.5*(dat.var['V'][1:,jy,:]+dat.var['V'][:-1,jy,:])
dpsx = (pspad[2:]-pspad[:-2])/(2*dmesh*np.cos(np.deg2rad(dat.lats[jy])))
dpsy = (dat.var['SP'][jy+1,:]-dat.var['SP'][jy-1,:])/(2*dmesh)
sig = (dat.var['PT'][:-1,jy,:]-dat.var['PT'][1:,jy,:]) \
     /(dat.var['P'][:-1,jy,:]-dat.var['P'][1:,jy,:])
tedges2 = 0.5 * (dat.var['PT'][:,jy,1:] + dat.var['PT'][:,jy,:-1])
wast = 86400*dat.var['W'][:-1,jy,:] * (dat.var['PT'][1:,jy,:]-dat.var['PT'][:-1,jy,:]) \
                             / (dat.var['P'][1:,jy,:]-dat.var['P'][:-1,jy,:])
tpad = np.pad(dat.var['PT'][:,jy,:],(1,1),'wrap')
dtx = (tpad[1:-1,2:]-tpad[1:-1,:-2])/(2*dmesh*np.cos(np.deg2rad(dat.lats[jy])))
dty = (dat.var['PT'][:,jy+1,:]-dat.var['PT'][:,jy-1,:])/(2*dmesh)
advt = dat.var['U'][:,jy,:] * dtx + dat.var['V'][:,jy,:] * dty
wast2 = wast + 0.5*(advt[1:,:]+advt[:-1,:])*86400

wast3 = wast2 - 86400 * bl[:,np.newaxis] * sig \
       *(U * dpsx[np.newaxis,:] + V * dpsx[np.newaxis,:])
# show wast3
ax = fig.add_subplot(223)       
im = ax.pcolormesh(xedges[ixr[0]-1:ixr[1]],tedges2[ktr[0]-1:ktr[1],ixr[0]-1:ixr[1]],
                wast3[ktr[0]:ktr[1],ixr[0]:ixr[1]],cmap=mymap,clim=(-3,5))
plt.clim(-3,5)
ax.contour(np.broadcast_to(dat.lons[ixr[0]:ixr[1]],(ktr[1]-ktr[0],ixr[1]-ixr[0])),
            tedges[ktr[0]:ktr[1],ixr[0]:ixr[1]],wast3[ktr[0]:ktr[1],ixr[0]:ixr[1]],[0,])
plt.ylim(330,400)
plt.xlim(30,120)
plt.tick_params(labelsize=fs)
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=fs)
cbar.ax.set_title('K/day',size=fs)
plt.xlabel('longitude',fontsize=fs)
plt.ylabel('potential temperature (K)',fontsize=fs)
plt.title(r'Vertical velocity as D$\theta$/Dt at 30N',fontsize=fs)

jy = 210
pspad = np.pad(dat.var['SP'][jy,:],(1,1),'wrap')
bl = 0.5*(dat.bl[1:]+dat.bl[:-1])
U = 0.5*(dat.var['U'][1:,jy,:]+dat.var['U'][:-1,jy,:])
V = 0.5*(dat.var['V'][1:,jy,:]+dat.var['V'][:-1,jy,:])
dpsx = (pspad[2:]-pspad[:-2])/(2*dmesh*np.cos(np.deg2rad(dat.lats[jy])))
dpsy = (dat.var['SP'][jy+1,:]-dat.var['SP'][jy-1,:])/(2*dmesh)
sig = (dat.var['PT'][:-1,jy,:]-dat.var['PT'][1:,jy,:]) \
     /(dat.var['P'][:-1,jy,:]-dat.var['P'][1:,jy,:])
tedges2 = 0.5 * (dat.var['PT'][:,jy,1:] + dat.var['PT'][:,jy,:-1])
wast = 86400*dat.var['W'][:-1,jy,:] * (dat.var['PT'][1:,jy,:]-dat.var['PT'][:-1,jy,:]) \
                             / (dat.var['P'][1:,jy,:]-dat.var['P'][:-1,jy,:])
tpad = np.pad(dat.var['PT'][:,jy,:],(1,1),'wrap')
dtx = (tpad[1:-1,2:]-tpad[1:-1,:-2])/(2*dmesh*np.cos(np.deg2rad(dat.lats[jy])))
dty = (dat.var['PT'][:,jy+1,:]-dat.var['PT'][:,jy-1,:])/(2*dmesh)
advt = dat.var['U'][:,jy,:] * dtx + dat.var['V'][:,jy,:] * dty
wast2 = wast + 0.5*(advt[1:,:]+advt[:-1,:])*86400
wast3 = wast2 - 86400 * bl[:,np.newaxis] * sig \
       *(U * dpsx[np.newaxis,:] + V * dpsx[np.newaxis,:])
# show wast3
ax = fig.add_subplot(224)       
im = ax.pcolormesh(xedges[ixr[0]-1:ixr[1]],tedges2[ktr[0]-1:ktr[1],ixr[0]-1:ixr[1]],
                wast3[ktr[0]:ktr[1],ixr[0]:ixr[1]],cmap=mymap,clim=(-3,5))
plt.clim(-3,5)
ax.contour(np.broadcast_to(dat.lons[ixr[0]:ixr[1]],(ktr[1]-ktr[0],ixr[1]-ixr[0])),
            tedges[ktr[0]:ktr[1],ixr[0]:ixr[1]],wast3[ktr[0]:ktr[1],ixr[0]:ixr[1]],[0,])
plt.ylim(330,400)
plt.xlim(30,120)
plt.tick_params(labelsize=fs)
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=fs)
cbar.ax.set_title('K/day',size=fs)
plt.xlabel('longitude',fontsize=fs)
plt.ylabel('potential temperature (K)',fontsize=fs)
plt.title(r'Vertical velocity as D$\theta$/Dt at 15N',fontsize=fs)

plt.show()

#%% Vertical velocity
fig = plt.figure()
ax = fig.add_subplot(111)
jy = 240
tedges2 = 0.5 * (dat.var['PT'][:,jy,1:] + dat.var['PT'][:,jy,:-1])
wast = 86400*dat.var['W'][:-1,jy,:] * (dat.var['PT'][1:,jy,:]-dat.var['PT'][:-1,jy,:]) \
                             / (dat.var['P'][1:,jy,:]-dat.var['P'][:-1,jy,:])
                          
ax = fig.add_subplot(223)                         
im = ax.pcolormesh(xedges[ixr[0]-1:ixr[1]],tedges2[ktr[0]-1:ktr[1],ixr[0]-1:ixr[1]],
                wast[ktr[0]:ktr[1],ixr[0]:ixr[1]],cmap=mymap,clim=(-3,5))
ax.contour(np.broadcast_to(dat.lons[ixr[0]:ixr[1]],(ktr[1]-ktr[0],ixr[1]-ixr[0])),
            tedges[ktr[0]:ktr[1],ixr[0]:ixr[1]],wast[ktr[0]:ktr[1],ixr[0]:ixr[1]],[0,])
plt.ylim(330,400)
plt.xlim(30,120)
plt.xlabel('longitude')
plt.ylabel('potential temperature (K)')
plt.title('Vertical velocity as heating rate DPT/Dt without advection (K/day)')
cbar = plt.colorbar(im)
plt.show()

# now apply the advective correction
# first step using only the derivative on model levels
# pad theta for the calculation of lon derivative
tpad = np.pad(dat.var['PT'][:,jy,:],(1,1),'wrap')
dtx = (tpad[1:-1,2:]-tpad[1:-1,:-2])/(2*dmesh*np.cos(np.deg2rad(dat.lats[jy])))
dty = (dat.var['PT'][:,jy+1,:]-dat.var['PT'][:,jy-1,:])/(2*dmesh)
advt = dat.var['U'][:,jy,:] * dtx + dat.var['V'][:,jy,:] * dty
wast2 = wast + 0.5*(advt[1:,:]+advt[:-1,:])*86400
# show wast2
ax = fig.add_subplot(224)
im = ax.pcolormesh(xedges[ixr[0]-1:ixr[1]],tedges2[ktr[0]-1:ktr[1],ixr[0]-1:ixr[1]],
                wast2[ktr[0]:ktr[1],ixr[0]:ixr[1]],cmap=mymap,clim=(-3,5))
ax.contour(np.broadcast_to(dat.lons[ixr[0]:ixr[1]],(ktr[1]-ktr[0],ixr[1]-ixr[0])),
            tedges[ktr[0]:ktr[1],ixr[0]:ixr[1]],wast2[ktr[0]:ktr[1],ixr[0]:ixr[1]],[0,])
plt.ylim(330,400)
plt.xlim(30,120)
plt.xlabel('longitude')
plt.ylabel('potential temperature (K)')
plt.title('Vertical velocity as heating rate DPT/Dt with advection (K/day)')
plt.colorbar(im)
plt.show()
# second step adding the correction due to variation of surface pressure
pspad = np.pad(dat.var['SP'][jy,:],(1,1),'wrap')
bl = 0.5*(dat.bl[1:]+dat.bl[:-1])
U = 0.5*(dat.var['U'][1:,jy,:]+dat.var['U'][:-1,jy,:])
V = 0.5*(dat.var['V'][1:,jy,:]+dat.var['V'][:-1,jy,:])
dpsx = (pspad[2:]-pspad[:-2])/(2*dmesh*np.cos(np.deg2rad(dat.lats[jy])))
dpsy = (dat.var['SP'][jy+1,:]-dat.var['SP'][jy-1,:])/(2*dmesh)
sig = (dat.var['PT'][:-1,jy,:]-dat.var['PT'][1:,jy,:]) \
     /(dat.var['P'][:-1,jy,:]-dat.var['P'][1:,jy,:])
wast3 = wast2 - 86400 * bl[:,np.newaxis] * sig \
       *(U * dpsx[np.newaxis,:] + V * dpsx[np.newaxis,:])
# show wast3
im = plt.pcolormesh(xedges[ixr[0]-1:ixr[1]],tedges2[ktr[0]-1:ktr[1],ixr[0]-1:ixr[1]],
                wast3[ktr[0]:ktr[1],ixr[0]:ixr[1]],cmap=mymap,clim=(-3,5))
plt.clim(-3,5)
plt.contour(np.broadcast_to(dat.lons[ixr[0]:ixr[1]],(ktr[1]-ktr[0],ixr[1]-ixr[0])),
            tedges[ktr[0]:ktr[1],ixr[0]:ixr[1]],wast3[ktr[0]:ktr[1],ixr[0]:ixr[1]],[0,])
plt.ylim(330,400)
plt.xlim(30,120)
plt.xlabel('longitude')
plt.ylabel('potential temperature (K)')
plt.title('Vertical velocity as heating rate DPT/Dt with full advection (K/day)')
plt.colorbar(im)
plt.show()

#%%
dat2 = dat.interpolPT(['W','ASPT','D','P','M'],[340,350,360,370,380,390,400],latrange=(0,50),lonrange=(-10,160))
#%%
dat3 = dat.interpolP(['W','U','V','D','ASPT','PT'],[23000,20000,17000,15000,13000,12000,11000,10000,9000,8000,7000],latrange=(0,50),lonrange=(-10,160))
#%%
dat2m = dat2.extract(lonrange=[67,93])
dat3m = dat3.extract(lonrange=[67,93])


