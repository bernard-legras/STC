#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Investigate the location of the high theta sources as seen in the diag data at high latitudes 

This script selects the sources with pressure less than 100 hPa. It counts them by 
latitude, longitude and level and calculate the mean pressure for each bin.

Then it plots the results after saving them in a file.

Created on Tue May  7 16:56:42 2019

@author: Bernard Legras
"""
import numpy as np
import pickle,gzip
from os.path import join
import matplotlib.pyplot as plt
from io107 import readpart107
import constants as cst
import socket
from numba import njit

streams = ["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"] 
#streams = ['Aug-21',]

# Find the main work dir according to the computer
if socket.gethostname() == 'gort':
    forw_dir =  '/home/legras/data/STC/STC-forw'
    init_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-initial'
elif 'ciclad' in socket.gethostname():
    forw_dir =  '/home/legras/STC/STC-forw'
    init_dir = '/data/legras/STC/STC-FORWBox-meanhigh-initial'
elif 'satie' in socket.gethostname():
    forw_dir =  '/home/legras/data/STC/STC-forw'
    init_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-initial'
elif 'Graphium' in socket.gethostname():
    forw_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-forw"
    init_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-FORWBox-meanhigh-initial" 
    
# area_pix: area of the pixel in the projected SAFNWC map (in degree^2)
area_pix = 0.1*0.1

# define vertical theta bins
minv = 322.5
maxv = 422.5
binv = 20
binvl = 5
vcent = np.arange(325,421,(maxv-minv)/binv)
#pedges = 100*np.arange(55,270,10)
#pcent = 100*np.arange(60,260,10)
# define horizontal grid
nlat = 50
nlon = 170
lats = np.arange(0.5,50.5)
lons = np.arange(-9.5,160.5)
xlon0 = -10
dlo = 1
ylat0 = 0
dla = 1
source_high = {}
source_high['count'] = np.zeros(shape=(binv,nlat,nlon),dtype=np.int32)
source_high['P'] = np.zeros(shape=(binv,nlat,nlon),dtype=np.float32)
source_high['lats'] = lats
source_high['lons'] = lons
source_high['vcent'] = vcent

@njit
def accumul(count,paccu,lev,jy,ix,pc):
    for i in range(len(lev)):
        count[lev[i],jy[i],ix[i]] += 1
        paccu[lev[i],jy[i],ix[i]] += pc[i]

for stream in streams:
    print(stream)
    # read the initial positions for thisi stream
    dir_in = join(init_dir,stream)
    pos0 = readpart107(0,dir_in,quiet=True)
     # extract the silviahigh selector
    ct = pos0['flag'] >> 24
    silviahigh = (ct == 9) | (ct == 13) | (ct == 8)
    # selection of the high clouds
    selec = silviahigh & (pos0['p']<27000)
    yc = pos0['y'][selec]
    xc = pos0['x'][selec]
    pc = pos0['p'][selec]
    tc = pos0['t'][selec]
    ir = pos0['ir_start'][selec]
    theta = tc * (cst.p0/pc)**cst.kappa
    ix = np.clip(np.floor((xc-xlon0)/dlo).astype(np.int),0,nlon-1)
    jy = np.clip(np.floor((yc-ylat0)/dla).astype(np.int),0,nlat-1)
    lev = np.clip(np.floor((theta-minv)/binvl).astype(np.int),0,binv-1)
    accumul(source_high['count'],source_high['P'],lev,jy,ix,pc)
    print('accumul',source_high['count'].max(),source_high['count'].min())
    # check the case with large theta at lats >40
    #selhigh = (yc<40) & (theta>380)
    # select the cases with low pressure (<100 hPa)
    # show their count as a function of time
    selhigh = (pc<=10000)
    print('nb high ',np.sum(selhigh),len(selhigh))
    hh,edges,_=plt.hist(ir[selhigh]/3600,bins=np.arange(-0.5,11*24+0.5))
    print(np.where(hh==0)[0])
    plt.show()
    
#%%
#=============================================================================
cc = source_high['count'].copy()
cc[cc==0] = 1
source_high['P'] /= cc
#=============================================================================

#=============================================================================
print ('backup')
with gzip.open(join(forw_dir,'source_high.pkl'),'wb') as f:
     pickle.dump(source_high,f)
#=============================================================================
#%%
# =============================================================================
fig = plt.figure(figsize=(15,12))
for lev in range(binv):
     plt.subplot(5,4,lev+1)
     im=plt.imshow(0.01*source_high['count'][lev,...],origin='lower',
                   extent=(-10,160,0,50),cmap='jet')    
     plt.colorbar(im,orientation='horizontal')
     plt.title(vcent[lev])
plt.show()
# =============================================================================
fig = plt.figure(figsize=(15,12))
for lev in range(binv):
     plt.subplot(5,4,lev+1)
     im=plt.imshow(0.01*source_high['count'][lev,...],origin='lower',
                   extent=(-10,160,0,50),cmap='jet')    
     plt.colorbar(im,orientation='horizontal')
     plt.title(vcent[lev])
plt.show()