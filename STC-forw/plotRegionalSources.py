# -*- coding: utf-8 -*-
"""
Created on 23 October 2019

Author: Bernard Legras
"""

from os.path import join
import pickle,gzip
import numpy as np
import matplotlib.pyplot as plt
from group import group

supertype = 'EID-FULL'
quiet = False
vert = 'theta'
target = 'global'
hmax = 1728

rea = {'EAD':'ERA5','EID-FULL':'EI','EAZ':'ERA5'}

out_dir = join('..','STC-FORWBox-meanhigh-OUT')

streams = ["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]

# read the mask
with gzip.open(join('..','mkSTCmask','MaskCartopy2-STCfine.pkl'),'rb') as f: 
    mm =pickle.load(f)
regcode = mm['regcode']
mask = mm['mask']
icx=icy=0.25
xlon0 = -10
ylat0 = 0
nlon = mm['nlons']
nlat = mm['nlats']

# target discretization (borrowed from transit.py)
minv = 322.5
maxv = 422.5
binv = 20
deltav = (maxv-minv)/binv
vcent = np.arange(minv+0.5*deltav,maxv,deltav)
vedge = np.arange(minv,maxv+0.01,deltav)

# source distribution in the same interval but with 200 bins
bins = 200 
deltas = (maxv-minv)/bins
scent = np.arange(minv+0.5*deltas,maxv,deltas)
sedge = np.arange(minv,maxv+0.01,deltas)

# read the first stream
run_type = supertype+'-Box-'+vert+'-'+target
regional_sav_stream = join(out_dir,'regionalN-save-stream-'+run_type+'-'+streams[0]+'-h'+str(hmax)+'.pkl')

with gzip.open(regional_sav_stream,'rb') as f:
    reg_src = pickle.load(f)
    
# read the other streams and accumulate in reg_src
for date in streams[1:]:
    regional_sav_stream = join(out_dir,'regionalN-save-stream-'+run_type+'-'+date+'-h'+str(hmax)+'.pkl')
    with gzip.open(regional_sav_stream,'rb') as f:
        buf = pickle.load(f)
        for reg in regcode:
            reg_src['count'][reg] += buf['count'][reg]
            reg_src['age'][reg] += buf['age'][reg]

for gr in group:
    reg_src['count'][gr] = np.zeros(shape=(binv,bins))
    reg_src['age'][gr] = np.zeros(shape=(binv,bins))
    for reg in group[gr]:
        reg_src['count'][gr] += reg_src['count'][reg]
        reg_src['age'][gr] += reg_src['age'][reg]
        
lzrh_file = join('..','..','MonthlyMeans','python','LZRH_Jul-Aug 2017.pkl')
with gzip.open(lzrh_file,'rb') as f:
    lzrh = pickle.load(f)
    
crossover_file = 'crossover.pkl'
with gzip.open(crossover_file,'rb') as f:
    cross = pickle.load(f)

# print the proportion of sources above the LZRH for the regions
# target level (11 for 380K)
levt = 11   
print(supertype,target,'target level',vcent[levt])
print('above LZRH diags')
print()
for reg in reg_src['count']:
    if np.isnan(lzrh[rea[supertype]][reg]['AS']):
        print('{:16} {:5.1f}%'.format(reg,100*np.sum(reg_src['count'][reg][levt,:])/np.sum(reg_src['count']['All'][levt,:])))
    else:
        lev0 = np.argmax(scent>lzrh[rea[supertype]][reg]['AS'])
        
        print('{:16} {:5.1f}% {:5.1f}%  {:5.1f}K'.format(reg,100*np.sum(reg_src['count'][reg][levt,:])/np.sum(reg_src['count']['All'][levt,:]),
                                  100*np.sum(reg_src['count'][reg][levt,lev0:])/np.sum(reg_src['count'][reg][levt,:]),lzrh[rea[supertype]][reg]['AS']))
del reg
print()
print('above crossover diags')
print()

# Do the same for the crossover and the groups
for gr in group.keys():
    if np.isnan(cross[target][supertype]['sh'][gr]):
        print('{:16} {:5.1f}%'.format(gr,100*np.sum(reg_src['count'][gr][levt,:])/np.sum(reg_src['count']['All'][levt,:])))
    else:
        lev0 = np.argmax(scent>cross[target][supertype]['sh'][gr])        
        print('{:16} {:5.1f}% {:5.1f}%  {:5.1f}K'.format(gr,100*np.sum(reg_src['count'][gr][levt,:])/np.sum(reg_src['count']['All'][levt,:]),
                                  100*np.sum(reg_src['count'][gr][levt,lev0:])/np.sum(reg_src['count'][gr][levt,:]),cross[target][supertype]['sh'][gr]))           