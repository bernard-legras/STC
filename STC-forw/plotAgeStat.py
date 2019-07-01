#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script provides statistics of the forward runs of the FORWBox-meanhigh family

It is based on the output of ageStat.py

It generates for both sh and mh hightypes 4  figures which are
(actually mh not plotted)

1) The 2d histogram of parcels as a function of age and potential temperature

2) The same histogram normalized per level to better see the vertical propagation

3) The mean age and the modal age for each level

4) The number of active parcels as a function of age

This is made for all of the 9 streams and for one supertype

EAZ, EAD, EIZ, EID are showing statistics of parcels in the FullAMA domain with the rule that a particle that exits
once is discarded.
EIZ-FULL, EID-FULL show the statistics of parcels in the global domain
EIZ-Return, EID-Return show the statistics of parcels in the FullAMA domain by applying a simple mask to the FULL runs, 
that is parcels that leave the domain are counted if they return.  

Created on Sat 23 Feb 2019

@author: Bernard Legras
"""
import numpy as np
#from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
import argparse
import socket
import os

parser = argparse.ArgumentParser()
parser.add_argument("-t","--type",type=str,choices=["EAD","EAZ","EIZ","EID","EIZ-Return","EID-Return","EIZ-FULL","EID-FULL"],help="type")
#parser.add_argument("-d","--date",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')

supertype = 'EIZ-FULL'
#hightypes = ['sh','mh']
hightypes = ['sh',]

step = 6
hmax = 1728
# 62 days
age_max = 1488

figsave = True

args = parser.parse_args()
if args.type is not None: supertype = args.type

dates = ["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]

# Find the main work dir according to the computer
if socket.gethostname() == 'gort':
    forw_dir =  '/home/legras/data/STC/STC-forw'
    out_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-OUT'
elif 'ciclad' in socket.gethostname():
    forw_dir =  '/home/legras/STC/STC-forw'
    out_dir = '/data/legras/STC/STC-FORWBox-meanhigh-OUT'
elif 'satie' in socket.gethostname():
    forw_dir =  '/home/legras/data/STC/STC-forw'
    out_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-OUT'
elif 'Graphium' in socket.gethostname():
    forw_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-forw"
    out_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-FORWBox-meanhigh-OUT" 

ns = int(hmax/step)
#fig = plt.figue(figsize=(13,13))
result = {}
for date in dates:
    print('')
    print('Processing '+date+' for '+supertype)
      
    # get the histograms and mean data calculated by checkmean.py from the part files
    file_out = os.path.join(out_dir,'ageStat-'+supertype+'-2017-'+date)
    print('reading',file_out )
    result[date] = {}
    result[date]['mh'] = {}
    result[date]['sh'] = {}
    with gzip.open(file_out,'rb') as f:
       [thets,result[date]['histog'],result[date]['edges'],result[date]['times6']] = pickle.load(f)
       for hightype in ['sh','mh']: 
           result[date][hightype]['nactiv'] = np.sum(result[date]['histog'][hightype],axis=1)

 #%% Plot of the age/thet histogram
for hightype in hightypes:
    fig = plt.figure(figsize=(11,10))
    fig.suptitle(hightype+' age theta histogram for '+supertype)
    n = 1
    for date in dates:
        plt.subplot(3,3,n)
        plt.imshow(np.log10(result[date]['histog'][hightype][0:248,:]).T,
                   extent=(0.25,62,275,600),origin='lower',aspect='auto',cmap='jet',clim=(-2,4))
        plt.ylim(300,450)
        plt.colorbar()
        plt.title(date)      
        if n in [1,4,7]: plt.ylabel('potential temperature (K)')
        if n in [7,8,9]: plt.xlabel('age (day)')
        n += 1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethethist-'+supertype))
    plt.show()
    
#%% Plot of the age/thet histogram, normalized for each level
for hightype in hightypes:
    fig = plt.figure(figsize=(11,10))
    fig.suptitle(hightype+' age theta normalized histogram for '+supertype)
    n = 1
    for date in dates:
        plt.subplot(3,3,n)
        hh = result[date]['histog'][hightype][0:248,:]
        ss = np.sum(hh,axis=0)
        hh = hh / ss[np.newaxis,:]
        plt.imshow(np.log10(hh.T),
                   extent=(0.25,62,275,600),origin='lower',aspect='auto',cmap='jet',clim=(-6,0))
        plt.ylim(300,450)
        plt.colorbar()
        plt.title(date)
        if n in [1,4,7]: plt.ylabel('potential temperature (K)')
        if n in [7,8,9]: plt.xlabel('age (day)')
        n += 1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agethetnormhist-'+supertype))
    plt.show()
    
#%% Mean age
ageaxis = np.arange(0.25,62.25,0.25) 
thetaxis = np.arange(275.5,700)
for hightype in hightypes:
    fig = plt.figure(figsize=(11,10))
    fig.suptitle(hightype+'  mean and modal age for '+supertype)
    n = 1
    for date in dates:
        plt.subplot(3,3,n)
        hh = result[date]['histog'][hightype][0:248,:]
        ss = np.sum(hh,axis=0)
        ss[ss==0]=1
        hh = hh / ss[np.newaxis,:]
        agemean = np.sum(hh*ageaxis[:,np.newaxis],axis=0)
        agemode = ageaxis[np.argmax(hh,axis=0)]
        plt.plot(agemean,thetaxis,'k',agemode,thetaxis,'r')
        plt.ylim(325,425)
        plt.title(date)
        if n in [1,4,7]: plt.ylabel('potential temperature (K)')
        if n in [7,8,9]: plt.xlabel('age (day)')
        n +=1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agemean-'+supertype))
    plt.show()

#%% Plot of the number of active parcels as a function of age
for hightype in hightypes:
    fig = plt.figure(figsize=(11,10))
    fig.suptitle(hightype+' nactiv for '+supertype+' as a function of age')
    n = 1
    for date in dates:
        plt.subplot(3,3,n)
        plt.semilogy(ageaxis,result[date][hightype]['nactiv'][0:248])
        plt.title(date)
        if n in [1,4,7]: plt.ylabel('nactiv')
        if n in [7,8,9]: plt.xlabel('age (day)')
        n += 1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs',hightype+'-agenactiv-'+supertype))
    plt.show()
