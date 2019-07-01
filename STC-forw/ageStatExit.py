#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script generates statistics on how the parcels escape from lower or lateral boundaries in the 
forward runs.The results are stored in the STC-FORWBox-meanhigh-OUT directory and need to be analysed 
by another program. The analysis must be done by reading the source location in the corresponding part_000
file (which is commun to all runs) and the end location which is stored in the output file of this script. 

Created on Mon 4 March 2019

@author: Bernard Legras
"""
import numpy as np
import pickle,gzip
from os.path import join
import argparse
import io107
#import constants as cst

#%%
parser = argparse.ArgumentParser()
#parser.add_argument("-y","--year",type=int,help="year")
#parser.add_argument("-m","--month",type=int,choices=1+np.arange(12),help="month")
parser.add_argument("-t","--type",choices=["EAD","EAZ","EIZ","EID","EIZ-Return","EID-Return","EIZ-FULL","EID-FULL"],help="type")
parser.add_argument("-d","--date",type=str,choices=["Jun-01","Jun-11","Jun-21","Jul-01","Jul-11","Jul-21","Aug-01","Aug-11","Aug-21"],help='run_date')
parser.add_argument("-i","--inc",type=int,help='step inc for processing')

# default values for the first start day
year = 2017
day = 1
supertype = 'EAZ'
step_inc = 6
date = 'Jun-01'

# pcut of the pressure
pcut = 45000

args = parser.parse_args()
if args.type is not None: supertype = args.type
if args.date is not None: date = args.date
if args.inc is not None: step_inc = args.inc

# ages max in days
age_max = 62
hmax = 24*(age_max + 10)
step_start = step_inc
max_injection_time = 11*24 + step_inc

#activate filtering if EID or EIZ
if supertype in ['EIZ','EID']:
    supertype1 = supertype + '-FULL'
    filtering = True
else:
    supertype1 = supertype
    filtering = False

# other parameters 
traj_dir = "/data/legras/flexout/STC/FORWBox-meanhigh"    
out_dir = "/data/legras/STC/STC-FORWBox-meanhigh-OUT"
file_out = join(out_dir,'ageStatExit-'+supertype+'-'+str(year)+'-'+date)

print("ageStatExit> process "+supertype+' '+date)
run_dir = join(traj_dir,'FORW-'+supertype1+'-'+str(year)+'-'+date)

# In order to have boxes every 6 h centered on 6, 12, ..., 1728
partStep = {}

# get the initial position to determine the weight
pos0 = io107.readpart107(0,run_dir,quiet=True)
y0 = pos0['y'].copy().astype(np.float32)
x0 = pos0['x'].copy().astype(np.float32)
p0 = pos0['p'].copy().astype(np.float32)
IDX_ORGN = pos0['idx_orgn']
numpart = pos0['numpart']
live = np.full(pos0['numpart'], fill_value = True, dtype=np.bool)
exited = {}
exited['exit'] = np.full(pos0['numpart'], fill_value = False, dtype=np.bool)
exited['bnd'] = np.full(pos0['numpart'], fill_value = 0, dtype=np.uint8)
exited['x'] = np.full(pos0['numpart'], fill_value = 0., dtype=np.float32)
exited['y'] = np.full(pos0['numpart'], fill_value = 0., dtype=np.float32)
exited['p'] = np.full(pos0['numpart'], fill_value = 0., dtype=np.float32)
exited['age'] = np.full(pos0['numpart'], fill_value = 0., dtype=np.float32)
#exited['ct'] =  (pos0['flag'] >> 24).astype(np.uint8)

startime = pos0['ir_start'].copy()/3600
idx_back0 = pos0['idx_back'].copy()

partStep[0] = {} 
partStep[0]['idx_back_init'] = idx_back0[startime <= step_start]
partStep[0]['idx_back'] = []

for step in range(step_start,hmax + step_inc ,step_inc):
    # Get rid of dictionary no longer used
    if step >= step_start + step_inc: del partStep[step-2*step_inc]
    # Read the new data
    partStep[step] = io107.readpart107(step,run_dir,quiet=True)
    # Set idx_back_init for next step
    if step <= max_injection_time:
        partStep[step]['idx_back_init'] = idx_back0[(startime > step) & (startime <= step + step_start)]
    else:
        partStep[step]['idx_back_init'] = None
    # Link the names as views
    data = partnow = partStep[step]
    partante = partStep[step-step_inc]
    #data['shfilt'] = silviahigh[selec]
    # filter applicable to EID-FULL and EIZ-FULL to get EID and EIZ cases
    # it may seem stupid to kill the parcel outside of the domain in the EI and EID case
    # to guess whee they have croddes later but this ensure the same method for ERA5 and ERA-I
    if filtering: 
        # kill the parcels which are out of the domain
        filterout = (data['x'] > 160) | (data['x']<-10) | (data['y']<0) | (data['y']>50)
        live[data['idx_back'][filterout]-IDX_ORGN] = False
        # selec the live parcels among the active ones 
        selec = live[data['idx_back']-IDX_ORGN]
        data['idx_back'] = data['idx_back'][selec]
        data['x'] = data['x'][selec]
        data['y'] = data['y'][selec]
        data['p'] = data['p'][selec]
        if data['idx_back_init'] is not None:
            data['idx_back_init'] = data['idx_back_init'][selec]

    # Make below a special case for FULL, in this case, we only search the exit through bottom boundary.
    # chek however there are no cases with exit throught the top

    # processing of deadborne parcels if any
    if partante['idx_back_init'] is not None:
        # get the index of the deadborne parcels
        borne = np.in1d(partante['idx_back_init'],partnow['idx_back'])
        idx_deadborne = partante['idx_back_init'][~borne]
        reborne = np.in1d(partnow['idx_back'],partante['idx_back_init'])
        print('# of new, borne, deadborne, reborne ',len(partante['idx_back_init']),np.sum(borne),len(idx_deadborne),np.sum(reborne))
        if len(idx_deadborne)>0:
            # set the parcels as exited and dead
            exited['exit'][idx_deadborne-IDX_ORGN] = True
            live[idx_deadborne-IDX_ORGN] = False
            if 'FULL' not in supertype:
                # find the closest lateral bndry (box 10W, 160E, 0, 50N ) in the order W, N, E, S
                proxyi = np.array([x0[idx_deadborne-IDX_ORGN]+10,50-y0[idx_deadborne-IDX_ORGN],
                                  160-x0[idx_deadborne-IDX_ORGN],y0[idx_deadborne-IDX_ORGN]]).argmin(axis=0)
                proxyv = np.array([x0[idx_deadborne-IDX_ORGN]+10,50-y0[idx_deadborne-IDX_ORGN],
                                  160-x0[idx_deadborne-IDX_ORGN],y0[idx_deadborne-IDX_ORGN]]).min(axis=0) 
                # difference with bottom boundary     
                proxp = pcut - p0[idx_deadborne-IDX_ORGN]
                # rescale and compare horizontal and vertical differences
                proxhv = np.array([proxp/1000,proxyv/3]).argmin(axis=0)
                # code the exitboundary (5 for bottom, 1 for W, 2 for N, 3 for E, 4 for S)
                exited['bnd'][idx_deadborne-IDX_ORGN][proxhv==0] = 5
                exited['bnd'][idx_deadborne-IDX_ORGN][proxhv==1] = 1+proxyi[proxhv==1]
            else:
                exited['bnd'][idx_deadborne-IDX_ORGN] = 5
            exited['x'][idx_deadborne-IDX_ORGN] = x0[idx_deadborne-IDX_ORGN]
            exited['y'][idx_deadborne-IDX_ORGN] = y0[idx_deadborne-IDX_ORGN]
            exited['p'][idx_deadborne-IDX_ORGN] = p0[idx_deadborne-IDX_ORGN]
            exited['age'][idx_deadborne-IDX_ORGN] = step_inc
    
    # processing the parcels which have exited between the two steps
    # only kept_a is to be calculated, kept_n is calculated for the purpose of checking    
    kept_a = np.in1d(partante['idx_back'],partnow['idx_back'],assume_unique=True)
    kept_n = np.in1d(partnow['idx_back'],partante['idx_back'],assume_unique=True)
    if partante['idx_back_init'] is not None:
        kept_n = kept_n & ~reborne      
    print('kept a, n ',len(kept_a),len(kept_n),kept_a.sum(),kept_n.sum(),(~kept_a).sum())
    
    if np.sum(~kept_a) >0:
        idx_lost = partante['idx_back'][~kept_a]
        # set the parcels as exited and dead
        exited['exit'][idx_lost-IDX_ORGN] = True
        live[idx_lost-IDX_ORGN] = False
        if 'FULL' not in supertype:
            # find the closest lateral bndry (box 10W, 160E, 0, 50N ) in the order W, N, E, S
            proxyi = np.array([partante['x'][~kept_a]+10,50-partante['y'][~kept_a],
                              160-partante['x'][~kept_a],partante['y'][~kept_a]]).argmin(axis=0)
            proxyv = np.array([partante['x'][~kept_a]+10,50-partante['y'][~kept_a],
                              160-partante['x'][~kept_a],partante['y'][~kept_a]]).min(axis=0)
            # difference with bottom boundary     
            proxp = pcut - partante['p'][~kept_a]
            # rescale and compare horizontal and vertical differences
            proxhv = np.array([proxp/1000,proxyv/3]).argmin(axis=0)
            # code the exitboundary (5 for bottom, 1 for W, 2 for N, 3 for E, 4 for S)
            exited['bnd'][idx_lost-IDX_ORGN][proxhv==0] = 5
            exited['bnd'][idx_lost-IDX_ORGN][proxhv==1] = 1+proxyi[proxhv==1]
        else:
            exited['bnd'][idx_deadborne-IDX_ORGN][partante['p'][~kept_a] > 10000] = 5
            exited['bnd'][idx_deadborne-IDX_ORGN][partante['p'][~kept_a] < 10000] = 6
        exited['x'][idx_lost-IDX_ORGN] = partante['x'][~kept_a]
        exited['y'][idx_lost-IDX_ORGN] = partante['y'][~kept_a]
        exited['p'][idx_lost-IDX_ORGN] = partante['p'][~kept_a]
        exited['age'][idx_lost-IDX_ORGN] = step - startime[idx_lost-IDX_ORGN]
  
with gzip.open(file_out,'wb') as f:
    pickle.dump(exited,f)
