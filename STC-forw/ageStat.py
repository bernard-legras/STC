#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script generate statistics of the distribution of parcels as a function of
the age with respect to convection. It is derived from checkmean.
The calculation is made by doing a two dimensional histogram in the age and the 
potential temperature along the trajectory.
It is done over the FullAMA region without separating the regions.
Weigthing is done by the factor (tau x pixel area/ Delta theta) factor
where tau is the interval between two images of the source clouds, 
pixel area is the pixel size in the PTOP file (used in prepforw5Box) 
and Delta theta is the vertical dicretization in theta.
tau is 1h and therefore does not appear explicitely since time unit is hour
(a factor 1/24 is introduced in plotting code to convert to day)
the pixel size is 0.1 x 0.1 degree^2 x cos(lat source)
Delta theta = 1K and therefore the factor does not appear

The wrong data on 30 August at 11h are filtered out.

The spurious sources at high latitude and high altitudes (theta > 360 and lat >40)
are filtered out.

Parameters:
--types: "EAD", "EAZ", "EIZ", "EID","EIZ-Return","EID-Return","EIZ-FULL","EID-FULL"
The first four are for trajectories remaining in the FullAMA domain. In EID-Return and
EIZ-Return, the trajectories are accounted when they are in the full domain but with a 
possibility of return. In EID-FULL ad EIZ-FULL, the statistics are performed in the
global domain.
--date: stream among "Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01",
"Aug-11","Aug-21"  
--inc: choice of the step increment

The output is contained in two dictionaries: 
    - thets contains basic information on the run
    min, max, mean (sh and mh) for all the steps
    - histog contains the mh and sh histograms
    Histograms are performed for theta bins between 275K and 700K (425 bins)
    and age bins of width step_inc and number hmax/step_inc (assumed to be an integer)
        
Stats are performed both for the meanhigh (mh) and silviahigh (sh) cloud tops. 

Output is stored in pkl file ageStat-xxx

Only one value of maximum age is used: 62 days = 1488 hours.
This generates hmax = 1728 hours because each stream is initialized over 10 days.

The checkmean script performs the same calculation but the histograms are made
as a function of time an not age. See the pograms in checkmean directory.

Created on Wed 20 Feb 2019

@author: Bernard Legras
"""
import numpy as np
import pickle,gzip
from os.path import join
import argparse
import io107
import constants as cst

#%%
parser = argparse.ArgumentParser()
#parser.add_argument("-y","--year",type=int,help="year")
#parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
parser.add_argument("-t","--type",choices=["EAT","EAD","EAZ","EIZ","EID","EIZ-Return","EID-Return","EIZ-FULL","EID-FULL"],help="type")
parser.add_argument("-d","--date",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')
parser.add_argument("-i","--inc",type=int,help='step inc for processing')
parser.add_argument("-c","--core",type=str,choices=["y","n"],help="core AMA (y) or not (n)")

# default values for the first start day
year = 2017
day = 1
supertype = 'EAD'
step_inc = 6
date = 'Jun-01'
core = False

# area_pix: area of the pixel in the projected SAFNWC map (in degree^2)
area_pix = 0.1*0.1

args = parser.parse_args()
if args.type is not None: supertype = args.type
if args.date is not None: date = args.date
if args.inc is not None: step_inc = args.inc
if args.core is not None:
    if args.core=='y': core=True
    else: core=False

# ages max in days
age_max = 62
hmax = 24*(age_max + 10)
step_start = step_inc

#activate filtering if EID or EIZ
if supertype in ['EIZ','EID']:
    supertype1 = supertype + '-FULL'
    filtering = True
    filtering_with_return = False
elif supertype in ['EIZ-Return','EID-Return']:
    supertype1 = supertype[:3] + '-FULL'
    filtering = False
    filtering_with_return = True
else:
    supertype1 = supertype
    filtering = False
    filtering_with_return = False

# other parameters 
traj_dir = "/data/legras/flexout/STC/FORWBox-meanhigh"    
out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"
file_out = join(out_dir,'ageStat-'+supertype+'-'+str(year)+'-'+date)
if core:
    file_out = join(out_dir,'ageStat-'+supertype+'-Core-'+str(year)+'-'+date)

print("ageStat> process "+supertype+' '+date)
run_dir = join(traj_dir,'FORW-'+supertype1+'-'+str(year)+'-'+date)

# define storage data
nbins = 425
edges = np.arange(275,701)
range_theta = [275,700]
times = np.arange(step_start,hmax+step_inc,step_inc)
nstep = int(hmax/step_inc)
# In order to have boxes every 6 h centered on 6, 12, ..., 1728
range_age = [3,1731]
histog = {}
histog['mh'] = np.zeros(shape=(nstep,nbins),dtype=np.float)
histog['sh'] = np.zeros(shape=(nstep,nbins),dtype=np.float)
thets = {}
thets['min_mh'] = np.empty(nstep,dtype=np.float32)
thets['max_mh'] = np.empty(nstep,dtype=np.float32)
thets['mean_mh'] = np.empty(nstep,dtype=np.float32)
thets['min_sh'] = np.empty(nstep,dtype=np.float32)
thets['max_sh'] = np.empty(nstep,dtype=np.float32)
thets['mean_sh'] = np.empty(nstep,dtype=np.float32)
nactiv = np.empty(nstep,dtype=np.int32)

# get the initial position to determine the weight
pos0 = io107.readpart107(0,run_dir,quiet=True)

y0 = pos0['y'].copy()
IDX_ORGN = pos0['idx_orgn']
numpart = pos0['numpart']
live = np.full(pos0['numpart'], fill_value = True, dtype='bool')
ct = pos0['flag'] >> 24
silviahigh = (ct == 9) | (ct == 13) | (ct == 8)
startime = pos0['ir_start'].copy()/3600
# filtering the 30 August at 11:00
if date == 'Aug-21':
    live[pos0['ir_start'] == 3600*227] = False
#  filtering high lat / high alt (spurious sources above 360K at lat>40N)
thet = pos0['t'] * (cst.p0/pos0['p'])**cst.kappa
live[(thet>360) & (pos0['y']>=40)] = False
# filtering out of the AMA core region if required
if core:
    live[(pos0['y']>40) | (pos0['y']<10) | (pos0['x']<20) | (pos0['x']>140)] = False
del thet
del pos0
del ct
thets['numpart'] = numpart
thets['bins'] = nbins
thets['range'] = range_theta
thets['step'] = step_inc
thets['hmax'] = hmax
thets['range_age'] = range_age
thets['nstep'] = nstep

l=0
for step in range(step_start,hmax + step_inc ,step_inc):
    data = io107.readpart107(step,run_dir,quiet=True)
    data['thet'] = data['t'] * (cst.p0/data['p'])**cst.kappa
    y0t = y0[data['idx_back']-IDX_ORGN]   
    ages = step - startime[data['idx_back']-IDX_ORGN]
    # no more need to filter the age
    #agefilt = ages <= age_max
    #shfilt = silviahigh[data['idx_back']-IDX_ORGN] & agefilt
    shfilt = silviahigh[data['idx_back']-IDX_ORGN]
    if filtering:
        # here we kill all parcels outside the FullAMA domain
        # designed for FULL runs with FullAMA target, that is supertype = EID or EIZ
        filterout = (data['x'] > 160) | (data['x']<-10) | (data['y']<0) | (data['y']>50)
        live[data['idx_back'][filterout]-IDX_ORGN] = False
        #selec = live[data['idx_back']-IDX_ORGN] & agefilt
        selec = live[data['idx_back']-IDX_ORGN]
        data['thet'] = data['thet'][selec]
        y0t = y0t[selec]
        shfilt = shfilt[selec]
        ages = ages[selec]
    elif filtering_with_return:
        # here we consider only the parcels in the FullAMA domain which can have visited the global domain
        # designed for the FULL runs with return target, that is supertype = EID-Return, EIZ-Return
        filterin = (data['x'] <= 160) & (data['x']>=-10) & (data['y']>=0) & (data['y']<=50) & live[data['idx_back']-IDX_ORGN]
        #data['thet'] = data['thet'][filterin & agefilt]
        #y0t = y0t[filterin & agefilt]
        #shfilt = shfilt[filterin & agefilt]
        data['thet'] = data['thet'][filterin]
        y0t = y0t[filterin]
        shfilt = shfilt[filterin]
        ages = ages[filterin]        
    else:
        # live filtering for other cases, that is either EAZ and EAD in the FullAMA domain or EIZ-FULL and EID-FULL
        # in the global domain
        selec = live[data['idx_back']-IDX_ORGN]
        data['thet'] = data['thet'][selec]
        y0t = y0t[selec]
        shfilt = shfilt[selec]
        ages = ages[selec]
        
    nactiv[l] = len(data['thet'])
    # try is used to face the case with no remaining parcels
    # stats for mh    
    try:
        print("thets mh",step,len(data['thet']),np.min(data['thet']),np.max(data['thet']),np.mean(data['thet']))
        histog['mh'] += np.histogram2d(ages,data['thet'],weights=area_pix * np.cos(np.deg2rad(y0t)),
               bins=[nstep,nbins],range=[range_age,range_theta])[0]
        print(np.sum(histog['mh']))
        thets['min_mh'][l] = np.min(data['thet'])
        thets['max_mh'][l] = np.max(data['thet'])
        thets['mean_mh'][l] = np.mean(data['thet'])
    except:
        print('uncompleted mh step',step)
    
    # stats for sh
    try:
        data['thet'] = data['thet'][shfilt]
        y0t = y0t[shfilt]
        ages = ages[shfilt]
        print("thets sh",step,len(data['thet']),np.min(data['thet']),np.max(data['thet']),np.mean(data['thet']))
        histog['sh'] += np.histogram2d(ages,data['thet'],weights=area_pix * np.cos(np.deg2rad(y0t)),
           bins=[nstep,nbins],range=[range_age,range_theta])[0]
        print(np.sum(histog['sh']))
        thets['min_sh'][l] = np.min(data['thet'])
        thets['max_sh'][l] = np.max(data['thet'])
        thets['mean_sh'][l] = np.mean(data['thet'])
    except: 
        print('uncompleted sh step',step)
    l += 1

histog['mh'] = histog['mh'].astype(np.float32)
histog['sh'] = histog['sh'].astype(np.float32)
with gzip.open(file_out,'wb') as f:
    pickle.dump([thets,histog,edges,times],f)
