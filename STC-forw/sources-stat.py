#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Establish stats on the source distributions
In the same style as in ageStat but from the part_000 files only
as part_000 is not read in ageStat
The output from this calculation can be used together with that of ageStat
to produce a plot of the histogram that includes age 0.
The calculations produces an histogram for the unions of all the streams
It uses as input the part_000 files contained in STC-FORWBox-meanhigh-initial.
The calculation is only performed for sh cloud tops.

Pruduces a file source_dist.pkl

Created on Sun May  5 21:38:53 2019

@author: Bernard Legras
"""
import numpy as np
import pickle,gzip
from os.path import join
import matplotlib.pyplot as plt
from io107 import readpart107
import constants as cst
import socket

streams = ["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]
#streams = ['Jul-21',]
# streams = ['Jun-01',]

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

# define vetical bins
nbins = 425
range_theta = [275,700]
source_dist = {}
source_dist['mh'] = np.zeros(nbins)
source_dist['sh'] = np.zeros(nbins)
source_dist['nbins'] = nbins
source_dist['range'] = range_theta
source_dist['thetac'] = np.arange(275.5,700.5)

#%%

for stream in streams:
    # read the initial positions for thisi stream
    dir_in = join(init_dir,stream)
    pos0 = readpart107(0,dir_in)
    y0 = pos0['y'].copy()
    theta0 = pos0['t'] * (cst.p0/pos0['p'])**cst.kappa
    IDX_ORGN = pos0['idx_orgn']
    numpart = pos0['numpart']
    print(stream, numpart)    
     # extract the silviahigh selector
    ct = pos0['flag'] >> 24
    silviahigh = (ct == 9) | (ct == 13) | (ct == 8)
    # filtering the 30 August at 11:00
    if stream == 'Aug-21':
        silviahigh = silviahigh & ~(pos0['ir_start'] == 3600*227)   
    source_dist['mh'] += np.histogram(theta0,bins=nbins,range=range_theta,weights=area_pix*np.cos(np.deg2rad(y0)))[0]
    source_dist['sh'] += np.histogram(theta0[silviahigh],bins=nbins,range=range_theta,weights=area_pix*np.cos(np.deg2rad(y0[silviahigh])))[0]
#%%    
with gzip.open(join(forw_dir,'source_dist.pkl'),'wb') as f:
    pickle.dump(source_dist,f)
#%%
with gzip.open(join(forw_dir,'source_dist.pkl'),'rb') as f:
    source_dist=pickle.load(f)
#%%
fs = 16
degree = 6372 * 2 * np.pi / 360
ff_s = (1/24) * degree**2 # conversion of time unit in da and degree into km
fig,ax = plt.subplots(figsize=(8,8))
ax.semilogx(ff_s*source_dist['sh'],source_dist['thetac'],'g',linewidth=8,alpha=0.5)
aa = source_dist['thetac'][np.argmax(source_dist['sh'])]
ax.plot([1.e4,1.e8],[aa,aa],'g')
ax.annotate(r'$\theta$ = 349.5 K',(1.50e4,351),fontsize=fs)
ax.tick_params(axis='x',labelcolor='g',labelsize=fs)
ax.set_xlim(1.e4,1.e8)
ax.set_ylim((320,440))
ax.set_xlabel(r'source distribution (day km$^2$ K$^{-1}$)',fontsize=fs,color='g')
ax.set_ylabel(r'potential temperature $\theta$ (K)',fontsize=fs)
ax.tick_params(labelsize=fs) 
plt.show()