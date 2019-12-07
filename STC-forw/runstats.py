#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script provides statistics of the forward runs of the FORWBox-meanhigh family

Should be applied with supertype = 'EAD', 'EAZ', 'EID-FULL', 'EIZ-FULL' 

It is based on the information contained in the traczilla log files.
The main technical problem is to read the sequence of these files, assumed to end as -a, -b, -c ,... 
for a given date, forgetting the missing ones. This is done imperfectly but still produces useful results.

It plots, as a function of time:
    - a figure showing the number of active parcels
    - a figure showing cpu usage
    - a figure showing zmin, zmax, zmean
    - a figure showing the 3 usefule components of npstop (type of exit from the domain)
    - a figure showing the 4 components of npsortie (exit side if lateral exit)

npstop
    2: (blue): out of the xylim domain (exit side in npsortie)
    3: (red): exit through the bottom boundary
    4: (cyan): exit through the upper boundary
    
npsortie:
    0: W (blue)
    1: N (red)
    2: E (black)
    3: S (cyan)

Created on Sat Jan 19 22:37:40 2019

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
import argparse
import socket
import os

parser = argparse.ArgumentParser()
parser.add_argument("-t","--type",type=str,choices=["EAD","EAZ","EIZ","EID","EIZ-FULL","EID-FULL"],help="type")
#parser.add_argument("-d","--date",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')

supertype = 'EAD'
hmax = 1728

step = 6

figsave = True

args = parser.parse_args()
if args.type is not None: supertype = args.type

dates = ["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"]

# Find the main work dir according to the computer
if socket.gethostname() == 'gort':
    main_work_dir  =  '/home/legras/data/STC/STC-work'
    forw_dir =  '/home/legras/data/STC/STC-forw'
    out_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-OUT'
elif 'ciclad' in socket.gethostname():
    main_work_dir = '/home/legras/flexpart/work/STC'
    forw_dir =  '/home/legras/STC/STC-forw'
    out_dir = '/data/legras/STC/STC-FORWBox-meanhigh-OUT'
elif 'satie' in socket.gethostname():
    main_work_dir = '/home/legras/data/STC/STC-work'
    forw_dir =  '/home/legras/STC/STC-forw'
    out_dir = '/home/legras/data/STC/STC-FORWBox-meanhigh-OUT'
elif socket.gethostname() == 'Graphium':
    main_work_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-work"
    forw_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-forw"
    out_dir = "C:\\cygwin64\\home\\berna\\data\\STC\\STC-FORWBox-meanhigh-OUT"
    
work_dir = os.path.join(main_work_dir,'FORWBox-meanhigh')
os.chdir(work_dir)
ll=os.listdir('.')

ns = int(hmax/step)
#fig = plt.figure(figsize=(13,13))
result = {}
for date in dates:
    print('')
    print('Processing '+date+' for '+supertype)
    # Generation of the main name of the file
    basicname = 'traczilla-FORW-'+supertype+'-2017-'+date 
    print('basicname',basicname)
    lls=[item for item in ll if basicname in item]
    if len(lls) ==0 : continue
    lls.sort()
    print(lls)
    time = np.arange(step,hmax+step,step)
    hours = np.zeros(ns,dtype=int)
    nactiv =  np.zeros(ns,dtype=int)
    npproc = np.zeros(ns,dtype=int)
    bounce = np.zeros(ns,dtype=int)
    cpu = np.zeros(ns,dtype=float)
    npstop = np.zeros((7,ns),dtype=int)
    npsortie = np.zeros((4,ns),dtype=int)
    zmin = np.zeros(ns,dtype=float)
    zmax = np.zeros(ns,dtype=float)
    zmean = np.zeros(ns,dtype=float)
    pmin = np.zeros(ns,dtype=float)
    pmax = np.zeros(ns,dtype=float)
    extension = False
    for file in lls:
        print('processing')
        if file != lls[0]:
            extension= True
            npsortie_base = npsortie[:,iso].copy()
        with open(file) as f:
            print('Opening',file)
            match = True
            cc = 0
            ff =  f.read().split('\n')
            l = -1
            while True:
                cc += 1
                if cc==300:
                    print('INFINITE LOOP DETECTED')
                    match = False
                    break
                try:
                    l +=1
                    line = ff[l]
                except:
                    print('EOF reached')
                    match = False
                    break
                if 'partout_stc' in line:
                    aa = line.split()
                    hour = np.int(aa[1].lstrip('part_'))
                    i = int((hour-step)/step)
                    cc = 0
                    #print('i',i)
                    # test the end of the run or the continuation
                    if i >= ns:
                        print('Reach beyond hmax')
                        match = False
                        break
                    hours[i] = hour
                    nactiv[i] = int(aa[2])
                    try:
                        l += 1
                        line = ff[l]
                    except:
                        print('unexpected EOF reached during incomplete sequence')
                        match = False
                        break
                    if 'checkpoint' in line:
                        l +=2 
                        line = ff[l]
                    if 'CONGRATULATIONS' in line: 
                        print('Normal end of run, incomplete sequence')
                        match = False
                        break
                    if 'npproc' not in line:
                        print('npproc not found',hour)
                        match = False
                        break
                    # decode npproc line
                    aa = line.split()
                    npproc[i] = int(aa[4])
                    bounce[i] = int(aa[5])
                    cpu[i] = float(aa[6].rstrip('s'))
                    # decode npstop line
                    aa = ff[l+1].split()
                    npstop[:,i] = np.array(aa[2:]).astype(np.int)
                    lc = l+1
                    # decode npsortie line for non FULL runs
                    if  'FULL' not in file:
                        lc += 1
                        aa = ff[lc].split()
                        if len(aa) == 5:
                            aa[4:6] = aa[4].split('E')
                        npsortie[:,i] = np.array([aa[2].rstrip(':W'),aa[3].rstrip(':N'),\
                                        aa[4].rstrip(':E'),aa[5].rstrip(':S')]).astype(np.int)                          
                        iso = i
                        if extension:
                            npsortie[:,i] += npsortie_base
                    # decode zmin, zmax, zmean
                    lc += 1
                    aa = ff[lc].split()
                    [zmin[i], zmax[i], zmean[i]] = [float(aa[1]), float(aa[2]), float(aa[3])]
                    lc += 1
                    # decode pmin, pmax
                    aa = ff[lc].split()
                    [pmin[i], pmax[i]] = [float(aa[1]), float(aa[2])]
                    l  = lc
                        
            if not match: continue              
         
    # store the data
    print('store the data')
    result[date] = {}
    result[date]['time'] = time
    result[date]['hours'] = hours
    result[date]['nactiv'] = nactiv 
    result[date]['npproc'] = npproc 
    result[date]['bounce'] = bounce
    # fix non meaningfull date in cpu
    cpu[cpu<0] = 0
    result[date]['cpu'] = cpu 
    result[date]['npstop'] = npstop 
    result[date]['npsortie'] = npsortie 
    result[date]['zmin'] = zmin
    result[date]['zmax'] = zmax
    result[date]['zmean'] = zmean
    result[date]['pmin'] = pmin
    result[date]['pmax'] = pmax

#%% plots of the diagnostics for the runs of this type
    
# Plot of the number of active parcels and npproc
fig = plt.figure(figsize=(9,9))
fig.suptitle('nactiv and npproc for '+supertype)
n = 1
for date in dates:
    plt.subplot(3,3,n)
    plt.semilogy(result[date]['time'],result[date]['nactiv'],'b',
             result[date]['time'],result[date]['npproc'],'r')
    plt.ylim(0,4e7)
    plt.title(date)
    n += 1
if figsave:
    plt.savefig(os.path.join(forw_dir,'figs','npproc-'+supertype))
plt.show()

#%% Plot of the cpu
fig = plt.figure(figsize=(9,9))
fig.suptitle('time step cpu (s) for '+supertype)
n = 1
for date in dates:
    plt.subplot(3,3,n)
    plt.plot(result[date]['time'],result[date]['cpu'],'b')
    #plt.ylim(-0.1e7,2.1e7)
    plt.title(date)
    n += 1
if figsave:
    plt.savefig(os.path.join(forw_dir,'figs','cpu-'+supertype))
plt.show()

 #%% Plot of the zmin, zmax, zmean
fig = plt.figure(figsize=(9,9))
fig.suptitle('zmin, max, zmean for '+supertype)
n = 1
for date in dates:
    plt.subplot(3,3,n)
    plt.plot(result[date]['time'],result[date]['zmin'],'b',
             result[date]['time'],result[date]['zmax'],'r',
             result[date]['time'],result[date]['zmean'],'k')
    if 'Z' in supertype:
        plt.ylim(-0.1,4)
    plt.title(date)
    n += 1
if figsave:
    plt.savefig(os.path.join(forw_dir,'figs','zminmax-'+supertype))
plt.show()

#%% Plot of the npstop
fig = plt.figure(figsize=(9,9))
fig.suptitle('npstop '+supertype)
n = 1
for date in dates:
    plt.subplot(3,3,n)
    plt.plot(result[date]['time'],result[date]['npstop'][0,:],'k',
             result[date]['time'],result[date]['npstop'][1,:],'g',
             result[date]['time'],result[date]['npstop'][2,:],'b',
             result[date]['time'],result[date]['npstop'][3,:],'r',
             result[date]['time'],result[date]['npstop'][4,:],'c',
             result[date]['time'],result[date]['npstop'][5,:],'m',
             result[date]['time'],result[date]['npstop'][6,:],'y')
    plt.title(date)
    n += 1
if figsave:
    plt.savefig(os.path.join(forw_dir,'figs','npstop-'+supertype))
plt.show()  

#%% Plot of the npsortie
if  'FULL'  not in supertype:
    fig = plt.figure(figsize=(9,9))
    fig.suptitle('npsortie '+supertype)
    n = 1
    for date in dates:
        plt.subplot(3,3,n)
        plt.plot(result[date]['time'],result[date]['npsortie'][0,:],'b',
                 result[date]['time'],result[date]['npsortie'][1,:],'r',
                 result[date]['time'],result[date]['npsortie'][2,:],'k',
                 result[date]['time'],result[date]['npsortie'][3,:],'c')
        plt.title(date)
        n += 1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs','npsortie-'+supertype))
    plt.show()
    
#%% Plot of the npsortie as diff
if  'FULL'  not in supertype:
    fig = plt.figure(figsize=(9,9))
    fig.suptitle('npsortie diff '+supertype)
    n = 1
    for date in dates:
        plt.subplot(3,3,n)
        diff = result[date]['npsortie'][:,1:] - result[date]['npsortie'][:,:-1]
        # truncate negatives values at discontinuities
        diff[diff<0] = 0
        plt.semilogy(result[date]['time'][1:],diff[0,:],'b',
                 result[date]['time'][1:],diff[1,:],'r',
                 result[date]['time'][1:],diff[2,:],'k',
                 result[date]['time'][1:],diff[3,:],'c')
        plt.title(date)
        n += 1
    if figsave:
        plt.savefig(os.path.join(forw_dir,'figs','npsortie-diff-'+supertype))
    plt.show()