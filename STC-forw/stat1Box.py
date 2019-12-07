# -*- coding: utf-8 -*-
"""

Produces statistics of the impact as a function of the vertical coordinate, collating data
from all the runs.
More precisely it sums the impact per 5K layer for mh and sh hightype for all the runs
and collect the result into a dictionary result which is then stored ito a pickle file.

This is intended for the comparison with calculations done with checkmean.py

Created on Fri Feb 15 2019

@author: Bernard
"""
import gzip, pickle
import transit as tt
#import argparse
import numpy as np
import socket
import os

hmax = 1728
vert = 'theta'
All = 'All'
water_path = True

types = ["EAZ","EIZ","EIZ-FULL","EAD","EID","EID-FULL"] 
dates = ["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]
#hightypes = ['mh','sh']

if socket.gethostname() == 'gort':
    out_dir = "/dkol/data/STC/STC-FORWBox-meanhigh-OUT"
    forw_dir =  '/home/legras/data/STC/STC-forw'
elif 'ciclad' in socket.gethostname():
    out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"
    forw_dir =  '/home/legras/data/STC/STC-forw'
elif socket.gethostname() == 'Graphium':
    out_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-FORWBox-meanhigh-OUT"
    forw_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-forw"
elif socket.gethostname() == 'satie':
    out_dir = "/limbo/data/STC/STC-FORWBox-meanhigh-OUT"
    forw_dir =  '/home/legras/data/STC/STC-forw'
else:
    'This program does not run on this computer'
    exit()
    
result = {}    

for type1 in types:
    
    result[type1] = {}
    if 'FULL' in type1:
        supertype = type1
        target = 'global'
    elif type1 in ['EIZ','EID']:
        supertype = type1 + '-FULL'
        target = 'FullAMA'
    else:
        supertype = type1
        target = 'FullAMA'

    run_type = supertype+'-Box-'+vert+'-'+target
 
    # Definition of the archive as a new transit class
    pile = tt.transit(water_path=water_path,vert=vert,target=target)
    
    # Accumulating the data from the multiple runs
    print("Accumulating for",type1)
    for date in dates:
        pile_save_stream = os.path.join(out_dir,'pile-save-stream-'+run_type+'-'+date+'-h'+str(hmax)+'.pkl')
        with gzip.open(pile_save_stream,'rb') as f:
            pile_arch = pickle.load(f)
            print(date,np.sum(pile_arch.transit['hist_t']),np.sum(pile_arch.transit['hist_t_vh']),np.sum(pile_arch.transit['hist_t_sh']))
        pile.merge(pile_arch)
        del pile_arch
    print(All,np.sum(pile.transit['hist_t']),np.sum(pile.transit['hist_t_vh']),np.sum(pile.transit['hist_t_sh']))
       
    # making averages
    #pile.complete()
    
    result[type1]['histsum_mh_t'] = np.sum(pile.transit['hist_t'],axis=(1,2))/5
    result[type1]['histsum_mh_s'] = np.sum(pile.transit['hist_s'],axis=(1,2))/5
    result[type1]['histsum_sh_t'] = np.sum(pile.transit['hist_t_sh'],axis=(1,2))/5
    result[type1]['histsum_sh_s'] = np.sum(pile.transit['hist_s_sh'],axis=(1,2))/5
    result[type1]['levsum'] = pile.vcent

with gzip.open(os.path.join(forw_dir,'stat1Box.pkl'),'wb') as f:
    pickle.dump(result,f)