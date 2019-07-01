#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Derived from sources-stat.py
In the same style as in ageStat but from the part_000 files only

In addition to sources-sta, performs an analysis of the distribution of sources
per region and groups of regions.

In the first part, generates a file source_dist_reg.pkl that contains these diagnostics.

In the second part, reads this file and prints diagnostic about the sources

The diagnostics are 
1) mean theta, theta value at the modal peak, percentage of sources from this 
region / total and / Asia
2) upward ratio for regional EAD crossover (proportion of sources above the local crossover)
   upward ratio for rmean EAD crossover (proportion of sources above the mean Asia crossover)
   upward ratio for regional EID crossover
   upward ratio for mean EID crossover

Created on Sun May  5 21:38:53 2019
Modified to add an analysis per region.

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
source_dist['nbins'] = nbins
source_dist['range'] = range_theta
source_dist['thetac'] = np.arange(275.5,700.5)

with gzip.open(join('..','mkSTCmask','MaskCartopy2-STCforwfine.pkl'),'rb') as f: 
    mm =pickle.load(f)
    
regcode = mm['regcode']
mask = mm['mask']
icx=icy=0.25
xlon0 = -10
ylat0 = 0
nlon = mm['nlons']
nlat = mm['nlats']

group = {}
group['Land'] = ['IndianSub','SouthChina','Pen','Pakistan','Bangladesh']
group['AsiaLand']  = group['Land']+['NorthChina','JapanKorea']
group['Seas'] = ['BoB','SCSPhi']
group['Ocean'] = group['Seas'] + ['IndianOcean','Indonesia','WestPacific','MidPacific']
group['Tibet'] = ['TibetanPlateau',]
group['AfricaArabia'] = ['CentralAfrica','GuineaGulf','NorthAfrica','RedSea','MiddleEast']
group['Asia'] = group['Ocean']+group['AsiaLand']+['NorthAsia','WestAsia','Caspian'] + ['TibetanPlateau',]
group['All'] = group['Asia'] + group['AfricaArabia'] + ['Europe','Atlantic','Mediterranea']

#%%
for reg in regcode.keys():
    source_dist[reg] = np.zeros(nbins)

for stream in streams:
    # read the initial positions for thisi stream
    dir_in = join(init_dir,stream)
    pos0 = readpart107(0,dir_in)
    y0 = pos0['y']
    x0 = pos0['x']
    idy = np.clip(np.floor((y0 - ylat0)/icy).astype(np.int),0,nlat-1)
    idx = np.clip(np.floor((x0 - xlon0)/icx).astype(np.int),0,nlon-1)
    reg_source = mask[idy,idx]
    theta0 = pos0['t'] * (cst.p0/pos0['p'])**cst.kappa
    IDX_ORGN = pos0['idx_orgn']
    numpart = pos0['numpart']
    #print(stream, numpart)
    #live = np.full(pos0['numpart'], fill_value = True, dtype='bool')
    # extract the silviahigh selector
    ct = pos0['flag'] >> 24
    silviahigh = (ct == 9) | (ct == 13) | (ct == 8)
    # filtering the 30 August at 11:00
    if stream == 'Aug-21':
        silviahigh = silviahigh & ~(pos0['ir_start'] == 3600*227) 
    
    for reg in regcode.keys():
        selec = (reg_source==regcode[reg]) & silviahigh
        source_dist[reg] += np.histogram(theta0[selec],bins=nbins,range=range_theta,weights=area_pix*np.cos(np.deg2rad(y0[selec])))[0]

for gr in group.keys():
    source_dist[gr] = np.zeros(nbins)
    for reg in group[gr]:
        source_dist[gr] += source_dist[reg]
        
#%%    
with gzip.open(join(forw_dir,'source_dist_reg.pkl'),'wb') as f:
    pickle.dump(source_dist,f)
    
#%% 
with gzip.open(join(forw_dir,'source_dist_reg.pkl'),'rb') as f:
    source_dist = pickle.load(f)

#%% stats sur les sources
gd_total = np.sum(source_dist['All'])
asia_total = np.sum(source_dist['Asia'])
for reg in list(regcode.keys()) + list(group.keys()):
    ss = np.sum(source_dist[reg])
    m1 = source_dist['thetac'][np.argmax(source_dist[reg])]
    m0 = np.sum(source_dist['thetac']*source_dist[reg])/np.sum(source_dist[reg])
    print('{:16} {:6.2f}K  {:6.2f}K {:5.1f}% {:5.1f}%'.format(reg,m0,m1,100*ss/gd_total,100*ss/asia_total))

#%%
fs = 16
ff_s = 1/24
fig,ax = plt.subplots(figsize=(8,8))
ax.semilogx(ff_s*source_dist['All'],source_dist['thetac'],'g',linewidth=8,alpha=0.5)
ax.tick_params(axis='x',labelcolor='g',labelsize=fs)
ax.set_xlim(1,4*10**4)
ax.set_ylim((320,440))
ax.set_xlabel(r'source distribution (day degree$^2$ K$^{-1}$)',fontsize=fs,color='g')
ax.set_ylabel(r'potential temperature $\theta$ (K)',fontsize=fs)
ax.tick_params(labelsize=fs) 
plt.show()

#%% Using the crossover and LZRH data, estimate the proportion of sources located
# above

# data from the LZRH and crossover estimations using EAD simulation sfor the crossover
# (MonthlyMeans/python/SRIP/MMBoxexNAll-Spe-CONF2019-2CICE.py  and diag-plot-2.py)
LZRH = {'TibetanPlateau':364.92,'AsiaLand':365.93,'Ocean':358.423,'Seas':354.92,'Asia':363.33}
CrossEAD = {'TibetanPlateau':364.16,'AsiaLand':365.41,'Ocean':362.5,'Seas':362.13,'Asia':363.90}
CrossEID = {'TibetanPlateau':363.1,'AsiaLand':361.9,'Ocean':359.1,'Seas':358.4,'Asia':361.9}
upward_ratio = {}
print('\nupward ratio for regional EAD crossover')
for reg in ['Asia','AsiaLand','Ocean','Seas','TibetanPlateau']:
    # proportion of sources above the local crossover
    selec = source_dist['thetac']>CrossEAD[reg]
    if reg == 'Asia':
        upward_ratio['Asia'] = np.sum(source_dist[reg][selec])/np.sum(source_dist[reg])
        print('{:16} upward ratio {:6.4f}'.format(reg,upward_ratio['Asia']))
    else:
        upward_ratio[reg] = np.sum(source_dist[reg][selec])/np.sum(source_dist[reg])
        factor = upward_ratio[reg] / upward_ratio['Asia']
        print('{:16} upward ratio {:6.4f} {:6.2f}'.format(reg,upward_ratio[reg],factor))
print('\nupward ratio for mean EAD crossover')
for reg in ['Asia','AsiaLand','Ocean','Seas','TibetanPlateau']:
    # proportion of sources above the mean Asia crossover
    selec = source_dist['thetac']>CrossEAD['Asia']
    if reg == 'Asia':
        upward_ratio['Asia'] = np.sum(source_dist[reg][selec])/np.sum(source_dist[reg])
        print('{:16} upward ratio {:6.4f}'.format(reg,upward_ratio['Asia']))
    else:
        upward_ratio[reg] = np.sum(source_dist[reg][selec])/np.sum(source_dist[reg])
        factor = upward_ratio[reg] / upward_ratio['Asia']
        print('{:16} upward ratio {:6.4f} {:6.2f}'.format(reg,upward_ratio[reg],factor))
print('\nupward ratio for regional EID crossover')
for reg in ['Asia','AsiaLand','Ocean','Seas','TibetanPlateau']:
    # proportion of sources above the local crossover
    selec = source_dist['thetac']>CrossEID[reg]
    if reg == 'Asia':
        upward_ratio['Asia'] = np.sum(source_dist[reg][selec])/np.sum(source_dist[reg])
        print('{:16} upward ratio {:6.4f}'.format(reg,upward_ratio['Asia']))
    else:
        upward_ratio[reg] = np.sum(source_dist[reg][selec])/np.sum(source_dist[reg])
        factor = upward_ratio[reg] / upward_ratio['Asia']
        print('{:16} upward ratio {:6.4f} {:6.2f}'.format(reg,upward_ratio[reg],factor))
print('\nupward ratio for mean EAD crossover')
for reg in ['Asia','AsiaLand','Ocean','Seas','TibetanPlateau']:
    # proportion of sources above the mean Asia crossover
    selec = source_dist['thetac']>CrossEID['Asia']
    if reg == 'Asia':
        upward_ratio['Asia'] = np.sum(source_dist[reg][selec])/np.sum(source_dist[reg])
        print('{:16} upward ratio {:6.4f}'.format(reg,upward_ratio['Asia']))
    else:
        upward_ratio[reg] = np.sum(source_dist[reg][selec])/np.sum(source_dist[reg])
        factor = upward_ratio[reg] / upward_ratio['Asia']
        print('{:16} upward ratio {:6.4f} {:6.2f}'.format(reg,upward_ratio[reg],factor))